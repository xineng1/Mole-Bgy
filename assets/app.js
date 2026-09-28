/* =============================================================
 * MolBench · app.js — 界面逻辑与可视化
 * ============================================================= */
(function () {
  'use strict';

  var C = window.MolCore, E = window.MolEnzymes, EX = window.MolExamples || [];
  function $(id) { return document.getElementById(id); }
  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }
  function num(v, d) { return C.fmt(v, d === undefined ? 2 : d); }

  /**
   * Markdown 表格单元格转义。
   * FASTA 头的名称里出现竖线非常常见（如 seq1|geneA|variant2），
   * 不转义会把一个单元格撑成多列，整张表错位。
   */
  function mdCell(s) {
    return String(s === null || s === undefined ? '' : s)
      .replace(/\|/g, '\\|')
      .replace(/\r?\n/g, ' ')
      .trim();
  }
  /** 分子量显示：≥10⁶ Da 改用 kDa，避免科学计数法不易读 */
  function fmtDa(v) {
    if (v === null || v === undefined || !isFinite(v)) return '—';
    if (Math.abs(v) >= 1e6) return num(v / 1000, 1) + ' kDa';
    return num(v, 1) + ' Da';
  }

  /* ---------------------------------------------------------
   * ORF 结果缓存
   * findORFs 是 O(6n) 的六框扫描，120 kb 序列单次约 45 ms。
   * 概览、位点图谱、环形图谱、报告、ORF 面板共 5 处都会用到，
   * 且切面板 / 改酶选择都会重绘 —— 不做缓存就会反复重算同一结果。
   * 用 seqId 做键（每次读入序列自增），比拿整条序列当键更省，也不会误命中。
   * --------------------------------------------------------- */
  var orfCache = new Map();

  function getORFs(minAa) {
    if (!state.seq || state.type === 'protein' || state.type === 'unknown') return [];
    var key = state.seqId + ':' + minAa + ':' + (state.circular ? 'c' : 'l');
    var hit = orfCache.get(key);
    if (hit) return hit;
    // 环状序列必须用环形扫描，否则会漏掉跨越首尾接缝的 ORF
    var r = state.circular
      ? C.findORFsCircular(state.seq, minAa, true)
      : C.findORFs(state.seq, minAa, true);
    if (orfCache.size >= 12) orfCache.clear();   // 防止长会话中无限增长
    orfCache.set(key, r);
    return r;
  }

  var state = {
    seq: '', name: '', type: 'dna', circular: false, features: [], seqId: 0,
    selectedEnzymes: []
  };

  var toastTimer = null;
  function toast(msg) {
    var t = $('toast');
    t.textContent = msg;
    t.classList.add('show');
    clearTimeout(toastTimer);
    toastTimer = setTimeout(function () { t.classList.remove('show'); }, 1800);
  }

  function copyText(text, what) {
    if (!text) { toast('没有可复制的内容'); return; }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(function () { toast((what || '内容') + '已复制'); },
        function () { fallbackCopy(text, what); });
    } else fallbackCopy(text, what);
  }
  function fallbackCopy(text, what) {
    var ta = document.createElement('textarea');
    ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
    document.body.appendChild(ta); ta.select();
    try { document.execCommand('copy'); toast((what || '内容') + '已复制'); }
    catch (e) { toast('复制失败，请手动选择'); }
    document.body.removeChild(ta);
  }

  /* =========================================================
   * Canvas 工具
   * ========================================================= */
  function setupCanvas(cv, cssHeight) {
    var dpr = window.devicePixelRatio || 1;
    var w = cv.clientWidth || cv.parentNode.clientWidth || 900;
    var h = cssHeight || cv.clientHeight || 300;
    cv.width = Math.round(w * dpr);
    cv.height = Math.round(h * dpr);
    cv.style.height = h + 'px';
    var ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx: ctx, w: w, h: h };
  }
  var COL = {
    grid: '#eef2f7', axis: '#cbd5e1', text: '#64748b', label: '#475569',
    gel: '#eef2f6', gelEdge: '#dbe3ec', band: '#1e3a8a', lane: '#ffffff',
    A: '#4f9d69', C: '#4c8dd9', G: '#e0a52b', T: '#d05a5a', U: '#a855f7', N: '#cbd5e1',
    line: '#2563eb', lineFill: 'rgba(37,99,235,.10)', ref: '#f59e0b'
  };

  /* =========================================================
   * 序列输入
   * ========================================================= */
  function readInput() {
    var parsed = C.parseFasta($('seqInput').value);
    state.seq = parsed.seq;
    state.name = parsed.name;
    state.type = parsed.type;
    state.circular = $('chkCircular').checked;
    state.features = parsed.features || [];
    state.seqId++;   // 序列一变，ORF 缓存键随之失效
    updateMeta();
    render(currentTab);
    histSchedule();
  }

  function updateMeta() {
    var s = state.seq;
    var warn = $('seqWarn');
    warn.classList.remove('show');
    if (!s) {
      $('seqMeta').innerHTML = '<span class="chip">尚未输入序列 — 可以点上方「载入示例序列」</span>';
      return;
    }
    var typeLabel = { dna: 'DNA', rna: 'RNA', protein: '蛋白质', unknown: '未知' }[state.type];
    var html = '';
    html += '<span class="chip">名称 <b>' + esc(state.name || '(未命名)') + '</b></span>';
    html += '<span class="chip">类型 <b>' + typeLabel + '</b></span>';
    html += '<span class="chip">长度 <b>' + s.length + '</b> ' + (state.type === 'protein' ? 'aa' : 'nt') + '</span>';
    if (state.type === 'dna' || state.type === 'rna') {
      html += '<span class="chip">GC <b>' + num(C.gcContent(s)) + '%</b></span>';
      html += '<span class="chip">拓扑 <b>' + (state.circular ? '环状' : '线性') + '</b></span>';
    }
    html += '<span class="chip">MW <b>' + fmtDa(getMw().primary) + '</b></span>';
    $('seqMeta').innerHTML = html;

    if (state.type === 'protein') {
      warn.textContent = '检测到的是蛋白质序列：序列操作、酶切、引物功能需要核酸序列。';
      warn.classList.add('show');
    } else if (state.type === 'unknown') {
      warn.textContent = '序列中含有无法识别的字符，请检查是否为有效的核酸或蛋白序列。';
      warn.classList.add('show');
    } else if (!C.isValidNucleic(s) && state.type !== 'protein') {
      warn.innerHTML = '序列包含简并碱基（R/Y/M/K/S/W/B/D/H/V/N）或非常规字符，酶切与翻译结果可能不完整。';
      warn.classList.add('show');
    }
  }

  /**
   * 分子量一次性算全并缓存。
   * 概览页会同时用到单链 / 双链 / RNA 三个值，加上 updateMeta 与统计磁贴，
   * 原先对同一条序列调用 5–6 次 mwNucleic —— 每次内部都要跑一遍 O(n) 的 baseCount，
   * 120 kb 下累计上百毫秒。现在只做一次计数。
   */
  var mwCache = { id: -1, val: null };

  function getMw() {
    if (mwCache.id === state.seqId && mwCache.val) return mwCache.val;
    var s = state.seq;
    var out;
    if (!s) {
      out = { ssDNA: 0, dsDNA: 0, ssRNA: 0, primary: 0 };
    } else if (state.type === 'protein') {
      out = { protein: C.mwProtein(s), primary: C.mwProtein(s) };
    } else {
      var c = C.baseCount(s);
      out = {
        ssDNA: C.mwNucleic(s, 'ssDNA', c),
        dsDNA: C.mwNucleic(s, 'dsDNA', c),
        ssRNA: C.mwNucleic(s, 'ssRNA', c)
      };
      out.primary = state.type === 'rna' ? out.ssRNA : out.dsDNA;
    }
    mwCache = { id: state.seqId, val: out };
    return out;
  }

  /* =========================================================
   * Tab 切换
   * ========================================================= */
  var currentTab = 'overview';
  Array.prototype.forEach.call(document.querySelectorAll('#tabs button'), function (btn) {
    btn.addEventListener('click', function () {
      Array.prototype.forEach.call(document.querySelectorAll('#tabs button'), function (b) { b.classList.remove('active'); });
      btn.classList.add('active');
      Array.prototype.forEach.call(document.querySelectorAll('.panel'), function (p) { p.classList.remove('active'); });
      var id = btn.getAttribute('data-panel');
      $('panel-' + id).classList.add('active');
      currentTab = id;
      render(id);
      window.scrollTo({ top: 0, behavior: 'smooth' });
    });
  });

  function render(tab) {
    if (tab === 'overview') renderOverview();
    else if (tab === 'seqops') renderSeqOps();
    else if (tab === 'digest') renderDigest();
    else if (tab === 'primer') renderPrimer();
  }

  /* =========================================================
   * 1. 概览
   * ========================================================= */
  function renderOverview() {
    var s = state.seq;
    if (!s) { $('ovStats').innerHTML = '<div class="stat"><div class="k">等待输入</div><div class="v">—</div></div>'; $('ovComp').innerHTML = ''; return; }

    var tiles = [];
    function tile(k, v, u) { tiles.push('<div class="stat"><div class="k">' + k + '</div><div class="v">' + v + (u ? '<span class="u">' + u + '</span>' : '') + '</div></div>'); }

    tile('长度', s.length, state.type === 'protein' ? ' aa' : ' nt');
    if (state.type !== 'protein') tile('GC 含量', num(C.gcContent(s)), '%');
    tile('分子量', fmtDa(getMw().primary), '');

    // 只跑一次 findORFs（原实现用 minAa=60 和 30 各跑一遍，长序列上白白翻倍）
    var orfsAll = getORFs(30);
    if (state.type !== 'protein') {
      var orfs = orfsAll.filter(function (o) { return o.aaLength >= 60; });
      tile('最长 ORF', orfs.length ? orfs[0].aaLength : 0, orfs.length ? ' aa' : '');
      tile('ORF 数量', orfs.length, ' 个');
      var mw = C.mwNucleic(s, 'dsDNA');
      tile('50 ng/µL 时', num(C.concNgUlToNm(50, mw), 1), ' nM');
    }
    var tmAll = (state.type === 'dna' && s.length >= 14 && s.length <= 60 && C.isValidNucleic(s))
      ? C.tmNearestNeighbor(s, {}) : NaN;
    tile('全序列 Tm', isNaN(tmAll) ? '—' : num(tmAll, 1), isNaN(tmAll) ? '' : ' °C');
    $('ovStats').innerHTML = tiles.join('');

    // 组成
    var cnt = C.baseCount(s);
    var total = s.length;
    var rows = '';
    Object.keys(cnt).forEach(function (k) {
      if (cnt[k] === 0) return;
      var pct = (cnt[k] / total) * 100;
      var color = COL[k] || '#94a3b8';
      var label = (state.type === 'protein') ? ('残基 ' + k) : (k === 'other' ? '其他字符' : (k === 'N' ? 'N / 简并' : k));
      rows += '<div style="display:flex;align-items:center;gap:10px;margin:5px 0">' +
        '<div style="width:92px;font-size:12.5px;color:var(--text-2)">' + esc(label) + '</div>' +
        '<div style="flex:1;background:var(--panel-2);border:1px solid var(--border);border-radius:5px;height:16px;overflow:hidden">' +
        '<div style="width:' + pct.toFixed(2) + '%;height:100%;background:' + color + '"></div></div>' +
        '<div style="width:130px;text-align:right;font-family:var(--mono);font-size:12px">' +
        cnt[k] + ' <span style="color:var(--text-3)">(' + num(pct, 1) + '%)</span></div></div>';
    });
    $('ovComp').innerHTML = rows;

    drawGC();
    drawBarcode();

    // 分子量
    $('mwSsDna').textContent = fmtDa(getMw().ssDNA);
    $('mwDsDna').textContent = fmtDa(getMw().dsDNA);
    $('mwRna').textContent = fmtDa(getMw().ssRNA);
    $('mwProt').textContent = orfsAll.length ? fmtDa(C.mwProtein(orfsAll[0].protein)) : '—';
    calcOD();
  }

  function drawGC() {
    var s = state.seq;
    var cv = $('gcCanvas');
    if (!s || state.type === 'protein') { var g0 = setupCanvas(cv, 220); g0.ctx.fillStyle = COL.text; g0.ctx.font = '13px system-ui'; g0.ctx.fillText('需要核酸序列', 12, 30); return; }
    var win = Math.max(5, parseInt($('gcWin').value, 10) || 100);
    var step = Math.max(1, parseInt($('gcStep').value, 10) || 10);
    var data = C.gcWindows(s, win, step);
    var g = setupCanvas(cv, 220);
    var ctx = g.ctx, W = g.w, H = g.h;
    var padL = 44, padR = 14, padT = 14, padB = 26;
    var pw = W - padL - padR, ph = H - padT - padB;

    // 网格
    ctx.strokeStyle = COL.grid; ctx.lineWidth = 1;
    ctx.font = '11px system-ui'; ctx.fillStyle = COL.text; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (var p = 0; p <= 100; p += 25) {
      var y = padT + ph - (p / 100) * ph;
      ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(W - padR, y); ctx.stroke();
      ctx.fillText(p + '%', padL - 6, y);
    }
    // 50% 参考线
    var y50 = padT + ph - 0.5 * ph;
    ctx.strokeStyle = COL.ref; ctx.setLineDash([4, 4]); ctx.beginPath();
    ctx.moveTo(padL, y50); ctx.lineTo(W - padR, y50); ctx.stroke(); ctx.setLineDash([]);

    if (!data.length) {
      ctx.fillStyle = COL.text; ctx.textAlign = 'center';
      ctx.fillText('序列长度小于窗口大小（' + win + ' bp）', W / 2, H / 2);
      $('gcCap').textContent = '';
      return;
    }
    var maxPos = data[data.length - 1].pos;
    var xs = function (pos) { return padL + ((pos - 1) / Math.max(1, maxPos - 1)) * pw; };
    var ys = function (gc) { return padT + ph - (gc / 100) * ph; };

    // 面积
    ctx.beginPath();
    ctx.moveTo(xs(data[0].pos), padT + ph);
    data.forEach(function (d) { ctx.lineTo(xs(d.pos), ys(d.gc)); });
    ctx.lineTo(xs(data[data.length - 1].pos), padT + ph);
    ctx.closePath();
    ctx.fillStyle = COL.lineFill; ctx.fill();

    ctx.beginPath();
    data.forEach(function (d, i) { i ? ctx.lineTo(xs(d.pos), ys(d.gc)) : ctx.moveTo(xs(d.pos), ys(d.gc)); });
    ctx.strokeStyle = COL.line; ctx.lineWidth = 1.8; ctx.stroke();

    // 坐标轴
    ctx.strokeStyle = COL.axis; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(padL, padT + ph); ctx.lineTo(W - padR, padT + ph); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(padL, padT); ctx.lineTo(padL, padT + ph); ctx.stroke();
    ctx.fillStyle = COL.text; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    var ticks = 5;
    for (var t = 0; t <= ticks; t++) {
      var pos = Math.round((maxPos * t) / ticks);
      ctx.fillText(pos + (t === ticks ? ' bp' : ''), xs(pos), padT + ph + 6);
    }
    var gcs = data.map(function (d) { return d.gc; });
    $('gcCap').textContent = '窗口 ' + win + ' bp / 步长 ' + step + ' bp，共 ' + data.length + ' 个窗口；' +
      'GC 范围 ' + num(Math.min.apply(null, gcs), 1) + '% ~ ' + num(Math.max.apply(null, gcs), 1) + '%，' +
      '均值 ' + num(gcs.reduce(function (a, b) { return a + b; }, 0) / gcs.length, 1) + '%';
  }

  function drawBarcode() {
    var s = state.seq;
    var g = setupCanvas($('barCanvas'), 60);
    var ctx = g.ctx, W = g.w, H = g.h;
    if (!s) return;
    var n = s.length;
    var cols = Math.min(n, Math.floor(W));
    var cw = W / cols;
    for (var i = 0; i < cols; i++) {
      var from = Math.floor((i * n) / cols), to = Math.max(from + 1, Math.floor(((i + 1) * n) / cols));
      var counts = {}, best = 'N', bestN = -1;
      for (var j = from; j < to && j < n; j++) {
        var ch = s.charAt(j);
        counts[ch] = (counts[ch] || 0) + 1;
        if (counts[ch] > bestN) { bestN = counts[ch]; best = ch; }
      }
      ctx.fillStyle = COL[best] || COL.N;
      ctx.fillRect(i * cw, 0, Math.ceil(cw), H);
    }
  }

  function calcOD() {
    var od = parseFloat($('odVal').value), dil = parseFloat($('odDil').value), kind = $('odKind').value;
    var s = state.seq;
    if (isNaN(od) || isNaN(dil)) { $('odConc').textContent = '—'; $('odNm').textContent = '—'; return; }
    var ugml = C.fromOD(od, dil, kind);
    $('odConc').textContent = num(ugml, 1) + ' µg/mL（= ' + num(ugml, 1) + ' ng/µL）';
    if (s && state.type !== 'protein') {
      var mw = kind === 'ssRNA' ? C.mwNucleic(s, 'ssRNA') : (kind === 'ssDNA' ? C.mwNucleic(s, 'ssDNA') : C.mwNucleic(s, 'dsDNA'));
      var nm = C.concNgUlToNm(ugml, mw);
      $('odNm').textContent = num(nm, 2) + ' nM = ' + num(nm / 1000, 3) + ' µM';
    } else $('odNm').textContent = '—';
  }

  /* =========================================================
   * 2. 序列操作
   * ========================================================= */
  function renderFeatures() {
    var card = $('featCard');
    var list = state.features || [];
    if (!list.length) { card.style.display = 'none'; return; }
    card.style.display = '';

    var html = '';
    list.forEach(function (f, i) {
      var label = f.qualifiers.gene || f.qualifiers.product || f.qualifiers.locus_tag || f.qualifiers.note || '—';
      if (f.qualifiers.product && f.qualifiers.gene && label === f.qualifiers.gene) {
        label = f.qualifiers.gene + ' — ' + f.qualifiers.product;
      }
      var tag = (f.type === 'CDS') ? '<span class="tag brand">CDS</span>'
        : (f.type === 'source' ? '<span class="tag">source</span>' : '<span class="tag">' + esc(f.type) + '</span>');
      html += '<tr><td>' + tag + '</td>' +
        '<td style="font-family:var(--mono);font-size:11.5px">' + esc(f.location) + '</td>' +
        '<td class="num">' + (f.start ? f.start + '–' + f.end : '—') + '</td>' +
        '<td>' + (f.strand === '-' ? '<span class="tag warn">−</span>' : '<span class="tag">+</span>') + '</td>' +
        '<td>' + esc(label) + '</td>' +
        '<td class="num">' + (f.length ? f.length + ' nt' : '—') + '</td>' +
        '<td>' + (f.start ? '<button class="small link" data-feat="' + i + '">查看</button>' : '') + '</td></tr>';
    });
    var tb = $('featTbl').querySelector('tbody');
    tb.innerHTML = html;
    Array.prototype.forEach.call(tb.querySelectorAll('button[data-feat]'), function (b) {
      b.addEventListener('click', function () {
        var f = (state.features || [])[parseInt(b.getAttribute('data-feat'), 10)];
        if (!f) return;
        var sub = C.extractFeature(state.seq, f);
        $('featTitle').textContent = '注释区段：' + f.type + ' ' + f.location +
          '（' + sub.length + ' nt' + (f.strand === '-' ? '，已按负链取反向互补' : '') + '）';
        $('featSeq').textContent = sub ? C.formatSeq(sub, 60) : '（该区段超出序列范围）';
        var note = [];
        if (f.qualifiers.translation) {
          note.push('GenBank 提供的翻译产物：' + f.qualifiers.translation.length + ' aa');
        }
        if (f.type === 'CDS' && sub.length % 3 === 0) {
          var mine = C.translate(sub);
          note.push('本地翻译：' + mine.replace(/\*$/, '') + (mine.charAt(mine.length - 1) === '*' ? '' : '（未以终止子结束）'));
        }
        $('featNote').textContent = note.join('；');
      });
    });
    $('featNote').textContent = '共 ' + list.length + ' 条注释。';
  }

  // 本页功能需要核酸序列；蛋白 / 未知类型时禁用按钮，避免产出全是 X 的垃圾翻译
  var SEQOPS_BTNS = ['btnRevComp', 'btnComp', 'btnRev', 'btnTranscribe', 'btnBackTranscribe', 'btnTranslate', 'btnFindOrf'];

  function renderSeqOps() {
    renderFeatures();
    var usable = !!state.seq && (state.type === 'dna' || state.type === 'rna');
    SEQOPS_BTNS.forEach(function (id) { $(id).disabled = !usable; });
    if (!usable) {
      var what = state.seq ? (state.type === 'protein' ? '蛋白质' : '未知类型') : '';
      $('opOut').textContent = state.seq
        ? '当前序列被识别为' + what + '，本页的反向互补 / 转录 / 翻译 / ORF 功能需要 DNA 或 RNA 序列。'
        : '点击上方按钮生成结果。';
      return;
    }
    if ($('opOut').textContent.indexOf('点击上方') === 0 ||
        $('opOut').textContent.indexOf('当前序列被识别为') === 0) {
      $('opOut').textContent = '点击上方按钮生成结果。';
    }
  }

  function showOp(text) {
    $('opOut').textContent = text || '（空）';
  }

  /* =========================================================
   * 3. 酶切
   * ========================================================= */
  function renderEnzymeList() {
    var kw = ($('enzSearch').value || '').trim().toUpperCase();
    var keepScroll = $('enzList').scrollTop;
    var html = '';
    E.list.forEach(function (en) {
      if (kw && en.name.toUpperCase().indexOf(kw) === -1 && en.siteClean.indexOf(kw) === -1) return;
      var checked = state.selectedEnzymes.indexOf(en.name) !== -1 ? ' checked' : '';
      html += '<label title="' + esc(en.site) + ' · ' + esc(en.note) + '">' +
        '<input type="checkbox" value="' + en.name + '"' + checked + '>' +
        en.name + ' <span class="site">' + esc(en.site) + '</span></label>';
    });
    $('enzList').innerHTML = html || '<span class="note">没有匹配的酶</span>';
    $('enzList').scrollTop = keepScroll;   // 重建列表后保持滚动位置

    // 只针对被点击的那一个酶做增删。
    // 早期实现是从 DOM 全量重建 selectedEnzymes，一旦处于搜索过滤状态，
    // 不在当前过滤结果里的已选酶（如搜索 "eco" 时的 HindIII）就会被静默丢弃。
    Array.prototype.forEach.call($('enzList').querySelectorAll('input'), function (inp) {
      inp.addEventListener('change', function () {
        var name = inp.value;
        var i = state.selectedEnzymes.indexOf(name);
        if (inp.checked) { if (i === -1) state.selectedEnzymes.push(name); }
        else if (i !== -1) state.selectedEnzymes.splice(i, 1);
        renderDigest();
      });
    });
  }

  /**
   * 注意：这里不再调用 renderEnzymeList()。
   * 酶列表有 51 个 checkbox，重建一次约 30 ms；而勾选某个酶、切换凝胶浓度、
   * 切到本面板都不会改变"列表内容本身"，重建纯属浪费。
   * 只有"选中集合被程序改动"（快捷选择）或搜索关键字变化时才需要重建，
   * 这些入口各自显式调用了 renderEnzymeList()。
   */
  function renderDigest() {
    var s = state.seq;
    var sel = state.selectedEnzymes.map(function (n) { return E.byName(n); }).filter(Boolean);

    if (!s || state.type === 'protein' || state.type === 'unknown') {
      $('siteTbl').querySelector('tbody').innerHTML = '<tr><td colspan="6" style="color:var(--text-3)">需要有效的核酸序列</td></tr>';
      $('fragTbl').querySelector('tbody').innerHTML = '';
      $('fragSummary').textContent = '';
      $('digestModeHint').textContent = '';
      var g = setupCanvas($('gelCanvas'), 440);
      g.ctx.fillStyle = COL.text; g.ctx.font = '13px system-ui';
      g.ctx.fillText('需要有效的核酸序列', 14, 30);
      mapHits = []; mapHover = -1; drawMap(null); drawPlasmidMap(null);
      return;
    }

    var d = C.digest(s, sel, state.circular);
    $('digestModeHint').textContent = (state.circular ? '当前按环状（质粒）计算' : '当前按线性计算') +
      ' · 序列 ' + s.length + ' bp';

    // 切点表
    var rows = '';
    var any = false;
    d.perEnzyme.forEach(function (pe) {
      pe.hits.forEach(function (h) {
        any = true;
        var cutOffset = pe.enzyme.cutOffset;
        var overhang = (cutOffset === Math.floor(pe.enzyme.length / 2) && pe.enzyme.length % 2 === 0) ? '平末端' : "5′ 突出";
        rows += '<tr><td><b>' + pe.enzyme.name + '</b></td>' +
          '<td style="font-family:var(--mono);font-size:12px">' + esc(pe.enzyme.site) + '</td>' +
          '<td class="num">' + h.start + '</td>' +
          '<td class="num">' + h.cut + ' / ' + (h.cut + 1) + '</td>' +
          '<td style="font-family:var(--mono);font-size:12px">' + esc(h.matched) + '</td>' +
          '<td><span class="tag">' + overhang + '</span>' +
          (h.wrap ? ' <span class="tag warn" title="该位点跨越序列首尾接缝，仅在环状模式下能识别">跨接缝</span>' : '') +
          '</td></tr>';
      });
    });
    if (!any) rows = '<tr><td colspan="6" style="color:var(--text-3)">' + (sel.length ? '所选酶在本序列中没有切点' : '请先在上方选择酶') + '</td></tr>';
    $('siteTbl').querySelector('tbody').innerHTML = rows;

    // 片段表
    var frows = '';
    if (sel.length) {
      var maxLen = d.fragments.length ? d.fragments[0] : 1;
      var rng = gelConf().range;
      d.fragments.forEach(function (f, i) {
        var intensity = Math.max(0.3, f / maxLen);
        var outNote = '';
        if (f > rng[1]) outNote = ' <span class="tag warn">大于分辨上限</span>';
        else if (f < rng[0]) outNote = ' <span class="tag warn">跑出胶外</span>';
        frows += '<tr><td class="num">' + (i + 1) + '</td>' +
          '<td class="num"><b>' + f + '</b>' + outNote + '</td>' +
          '<td class="num">' + num((f / s.length) * 100, 1) + '%</td>' +
          '<td class="num">' + num(migrationMm(f), 1) + ' mm</td>' +
          '<td><div style="width:120px;height:9px;background:var(--panel-2);border:1px solid var(--border);border-radius:3px">' +
          '<div style="width:' + (intensity * 100).toFixed(0) + '%;height:100%;background:' + COL.band + ';opacity:' + intensity.toFixed(2) + '"></div></div></td></tr>';
      });
    } else {
      frows = '<tr><td colspan="5" style="color:var(--text-3)">未选择酶，无切割</td></tr>';
    }
    $('fragTbl').querySelector('tbody').innerHTML = frows;
    $('fragSummary').innerHTML = sel.length
      ? '共 <b>' + d.fragments.length + '</b> 个片段，合计 <b>' + d.fragments.reduce(function (a, b) { return a + b; }, 0) + '</b> bp（原始 ' + s.length + ' bp）。' +
        '迁移距离为 1% 琼脂糖、约 90 V 恒压电泳的经验估算，仅供判断片段能否分开。'
      : '';

    drawGel(d, sel);
    drawMap(d);
    drawPlasmidMap(d);
  }

  /**
   * 各浓度琼脂糖的参数集中在一处，避免"分辨范围"与"迁移系数"两张表各写各的而失同步。
   * range = 有效分辨范围 (bp)；coef = 经验迁移系数，d(mm) = a − b·log₁₀(bp)
   */
  var GEL = {
    '0.7': { range: [400, 22000], coef: [70, 14.5] },
    '1':   { range: [200, 12000], coef: [60, 13.3] },
    '1.5': { range: [100, 8000],  coef: [50, 11.5] },
    '2':   { range: [80, 5000],   coef: [42, 10.0] }
  };
  function gelConf() {
    return GEL[String(parseFloat($('gelPct').value) || 1)] || GEL['1'];
  }

  function migrationMm(bp) {
    var coef = gelConf().coef;
    var d = coef[0] - coef[1] * Math.log10(bp);
    return Math.max(1.5, Math.min(coef[0], d));
  }

  var LADDER_1KB = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250];

  function drawGel(d, sel) {
    var cv = $('gelCanvas');
    var g = setupCanvas(cv, 440);
    var ctx = g.ctx, W = g.w, H = g.h;

    var perEnzyme = $('gelPerEnzyme').checked;
    var showLabel = $('gelLabel').checked;
    var pct = parseFloat($('gelPct').value) || 1;
    var range = gelConf().range;

    var lanes = [{ name: 'Ladder', frags: LADDER_1KB.slice(), isLadder: true }];
    if (state.seq) {
      if (sel.length === 0) {
        // 未酶切的环状质粒在胶上呈三种构象：超螺旋最快、线性居中、开环（切口）最慢。
        // 表观迁移长度按经验系数折算，从而落在不同的凝胶位置。
        if (state.circular) {
          lanes.push({
            name: '未酶切（质粒）',
            bands: [
              { bp: state.seq.length, apparent: state.seq.length * 0.70, label: 'SC' },
              { bp: state.seq.length, apparent: state.seq.length * 1.00, label: 'LIN' },
              { bp: state.seq.length, apparent: state.seq.length * 1.45, label: 'OC' }
            ]
          });
        } else {
          lanes.push({ name: '未酶切', frags: [state.seq.length] });
        }
      } else {
        lanes.push({ name: sel.map(function (e) { return e.name; }).join('+') || '—', frags: d.fragments.slice(), isMain: true });
        if (perEnzyme && sel.length > 1) {
          // 复用 renderDigest 已经算好的切点，只重算片段长度（O(k log k)），
          // 不再对每个酶把整条序列重新跑一遍正则
          d.perEnzyme.forEach(function (pe) {
            var cuts = pe.hits.map(function (h) { return h.cut; });
            lanes.push({
              name: pe.enzyme.name,
              frags: C.fragmentsFromCuts(state.seq.length, cuts, state.circular)
            });
          });
        }
      }
    }

    var padT = 34, padB = 46, padL = 52, padR = 18;
    var laneH = H - padT - padB;
    var gap = 8;
    var MIN_LANE_W = 18;

    // 按画布实际宽度决定最多画几条泳道。
    // 旧实现固定 14 条且把 laneW 硬 clamp 到 18，窄屏时总宽会超出画布右边界。
    var maxLanes = Math.max(2, Math.floor((W - padL - padR + gap) / (MIN_LANE_W + gap)));
    var truncated = lanes.length - maxLanes;
    if (truncated > 0) lanes = lanes.slice(0, maxLanes);
    var laneW = (W - padL - padR - gap * (lanes.length - 1)) / lanes.length;

    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);
    ctx.fillStyle = COL.gel;
    ctx.strokeStyle = COL.gelEdge; ctx.lineWidth = 1;
    ctx.fillRect(padL, padT, W - padL - padR, laneH);
    ctx.strokeRect(padL, padT, W - padL - padR, laneH);

    var logMin = Math.log10(range[0]), logMax = Math.log10(range[1]);
    function yOf(bp) {
      if (bp <= 0) return padT;
      var t = (logMax - Math.log10(bp)) / (logMax - logMin);
      t = Math.max(0, Math.min(1, t));
      return padT + t * laneH;
    }

    // Ladder 刻度
    ctx.font = '10px system-ui'; ctx.fillStyle = COL.text; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    LADDER_1KB.forEach(function (bp) {
      var y = yOf(bp);
      if (y < padT - 1 || y > padT + laneH + 1) return;
      ctx.strokeStyle = '#cbd5e1'; ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(padL - 4, y); ctx.lineTo(padL, y); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillText(bp >= 1000 ? (bp / 1000) + ' kb' : bp + ' bp', padL - 7, y);
    });

    lanes.forEach(function (lane, li) {
      var x0 = padL + li * (laneW + gap);
      // 加样孔
      ctx.fillStyle = '#e2e8f0';
      ctx.fillRect(x0 + laneW * 0.28, padT - 9, laneW * 0.44, 6);
      // 泳道底色
      ctx.fillStyle = 'rgba(255,255,255,.55)';
      ctx.fillRect(x0, padT, laneW, laneH);

      // 统一成 {bp, apparent, label}：普通片段 apparent === bp，构象带则按表观大小迁移
      var bands = lane.bands || lane.frags.map(function (bp) { return { bp: bp, apparent: bp }; });
      var maxF = Math.max.apply(null, bands.map(function (b) { return b.apparent; }).concat([1]));
      bands.forEach(function (band) {
        var bp = band.bp, app = band.apparent || bp, lbl = band.label;
        var y = yOf(app);
        if (app > range[1]) {
          // 过大：留在加样孔附近，画一条贴顶的淡带
          ctx.fillStyle = 'rgba(30,58,138,.45)';
          ctx.fillRect(x0 + 2, padT + 1.5, laneW - 4, 2.5);
          return;
        }
        if (app < range[0]) {
          // 过小：已跑出胶外，画一条贴底的淡带提示
          ctx.fillStyle = 'rgba(30,58,138,.35)';
          ctx.fillRect(x0 + 2, padT + laneH - 4, laneW - 4, 2.5);
          return;
        }
        var inten = lane.isLadder ? 0.75 : Math.max(0.28, Math.min(1, app / maxF));
        var bandH = lane.isLadder ? 3.2 : Math.max(3, Math.min(7, 3 + Math.log10(app)));
        ctx.fillStyle = 'rgba(30,58,138,' + inten.toFixed(2) + ')';
        ctx.fillRect(x0 + 2, y - bandH / 2, laneW - 4, bandH);
        if (showLabel && (lane.isLadder || lane.isMain || lanes.length <= 6)) {
          ctx.fillStyle = lbl ? '#7c3aed' : COL.label;
          ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
          ctx.font = (lbl ? 'bold ' : '') + '9.5px system-ui';
          var text = lbl ? lbl : (bp >= 1000 ? num(bp / 1000, 1) + 'k' : String(bp));
          ctx.fillText(text, x0 + laneW / 2, y - bandH / 2 - 1);
          if (lbl) {
            ctx.fillStyle = COL.text;
            ctx.font = '9px system-ui';
            ctx.textBaseline = 'top';
            ctx.fillText(bp >= 1000 ? num(bp / 1000, 1) + 'k' : String(bp), x0 + laneW / 2, y + bandH / 2 + 1);
          }
        }
      });

      // 泳道名（竖排旋转，避免拥挤）
      ctx.save();
      ctx.translate(x0 + laneW / 2, padT + laneH + 8);
      ctx.fillStyle = COL.label; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.font = (lane.isMain ? 'bold ' : '') + '10.5px system-ui';
      ctx.fillText(lane.name.length > 22 ? lane.name.slice(0, 21) + '…' : lane.name, 0, 0);
      ctx.restore();
    });

    // 方向标注
    ctx.fillStyle = COL.text; ctx.font = '10.5px system-ui'; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText('⊕ 电极方向：小片段在前（下方）', padL, 8);
    ctx.textAlign = 'right';
    ctx.fillText(pct + '% 琼脂糖 · 分辨范围约 ' + range[0] + '–' + range[1] + ' bp', W - padR, 8);
    if (truncated > 0) {
      ctx.fillStyle = '#b45309';
      ctx.textAlign = 'left';
      ctx.fillText('※ 泳道过多，仅显示前 ' + maxLanes + ' 条（另有 ' + truncated + ' 条未画出）', padL + 210, 8);
    }
  }

  /* =========================================================
   * 酶切位点线性图谱
   * ========================================================= */
  var MAP_COLORS = ['#2563eb', '#dc2626', '#059669', '#d97706', '#7c3aed',
                    '#0891b2', '#db2777', '#65a30d', '#ea580c', '#4f46e5'];
  var mapHits = [];      // 切点屏幕坐标，供悬停命中
  var mapHover = -1;
  var mapData = null;    // 最近一次绘制所用的数据，悬停重绘时复用

  function drawMap(d) {
    mapData = d;
    var s = state.seq;
    var cv = $('mapCanvas');
    if (!s || !d) { setupCanvas(cv, 130); return; }

    var enzymes = d.perEnzyme.filter(function (pe) { return pe.hits.length; });
    // 图例每 4 个一行，画布高度随行数增长
    var legendRows = Math.max(1, Math.ceil(enzymes.length / 4));
    var H = 96 + legendRows * 17;
    var g = setupCanvas(cv, H);
    var ctx = g.ctx, W = g.w;
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);

    var padL = 46, padR = 18;
    var axisY = 40, axisH = 16;
    var pw = W - padL - padR;
    var n = s.length;

    // 序列轴
    ctx.fillStyle = '#eef2f7';
    ctx.strokeStyle = COL.gelEdge; ctx.lineWidth = 1;
    roundRect(ctx, padL, axisY, pw, axisH, 4); ctx.fill(); ctx.stroke();

    // 刻度
    ctx.font = '10px system-ui'; ctx.fillStyle = COL.text;
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (var t = 0; t <= 4; t++) {
      var pos = Math.round((n * t) / 4);
      var x = padL + (pos / n) * pw;
      ctx.strokeStyle = '#cbd5e1';
      ctx.beginPath(); ctx.moveTo(x, axisY + axisH); ctx.lineTo(x, axisY + axisH + 4); ctx.stroke();
      ctx.fillText(pos >= 1000 ? (pos / 1000).toFixed(1) + 'k' : String(pos), x, axisY + axisH + 6);
    }
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillStyle = COL.text;
    ctx.fillText('序列', padL - 8, axisY + axisH / 2);

    // ORF 轨道（正链在上、负链在下）
    // ORF 查找对超长序列开销大，超过 200 kb 时跳过该轨道
    var orfs = (n <= 200000) ? getORFs(60).slice(0, 8) : [];
    var orfTop = axisY + axisH + 18;
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.font = '10px system-ui';
    ctx.fillStyle = COL.text; ctx.fillText('ORF', padL - 8, orfTop + 8);
    ctx.fillStyle = '#f1f5f9'; roundRect(ctx, padL, orfTop, pw, 16, 3); ctx.fill();
    orfs.forEach(function (o) {
      var x1 = padL + ((o.ntStart - 1) / n) * pw;
      var x2 = padL + (o.ntEnd / n) * pw;
      var y = o.strand === '+' ? orfTop + 2 : orfTop + 9;
      ctx.fillStyle = o.strand === '+' ? 'rgba(13,148,136,.75)' : 'rgba(124,58,237,.65)';
      ctx.fillRect(x1, y, Math.max(1.5, x2 - x1), 6);
    });

    // 切点竖线
    mapHits = [];
    enzymes.forEach(function (pe, ei) {
      var color = MAP_COLORS[ei % MAP_COLORS.length];
      pe.hits.forEach(function (h) {
        var x = padL + (h.cut / n) * pw;
        mapHits.push({ x: x, cut: h.cut, start: h.start, enzyme: pe.enzyme, matched: h.matched, wrap: h.wrap, color: color });
      });
    });
    mapHits.sort(function (a, b) { return a.x - b.x; });
    mapHits.forEach(function (hit, i) {
      var isHover = (i === mapHover);
      ctx.strokeStyle = hit.color;
      ctx.lineWidth = isHover ? 3.5 : 2;
      ctx.beginPath();
      ctx.moveTo(hit.x, axisY - 8);
      ctx.lineTo(hit.x, axisY + axisH + 3);
      ctx.stroke();
      if (isHover) {
        ctx.fillStyle = hit.color;
        ctx.beginPath(); ctx.arc(hit.x, axisY - 11, 3, 0, Math.PI * 2); ctx.fill();
      }
    });

    // 图例（与环形图谱共用同一套排版）
    drawLegend(ctx,
      enzymes.map(function (pe, ei) {
        return { color: MAP_COLORS[ei % MAP_COLORS.length], text: pe.enzyme.name + ' ×' + pe.hits.length };
      }),
      padL, orfTop + 26, pw);

    // 悬停提示
    if (mapHover >= 0 && mapHits[mapHover]) {
      var hit = mapHits[mapHover];
      var txt = hit.enzyme.name + '  ' + hit.enzyme.site + '  切点 ' + hit.cut + '/' + (hit.cut + 1) +
        '  识别起点 ' + hit.start + (hit.wrap ? '  [跨接缝]' : '');
      ctx.font = '11px system-ui';
      var tw = ctx.measureText(txt).width + 16;
      var bx = Math.min(Math.max(4, hit.x - tw / 2), W - tw - 4);
      var by = 4;
      ctx.fillStyle = 'rgba(31,41,55,.94)';
      roundRect(ctx, bx, by, tw, 22, 5); ctx.fill();
      ctx.fillStyle = '#fff'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText(txt, bx + 8, by + 11);
    }

    $('mapCap').textContent = enzymes.length
      ? '共 ' + mapHits.length + ' 个切点。ORF 轨道中上排为正链、下排为负链，仅显示最长的 8 个 ORF（≥60 aa）。'
      : '所选酶在当前序列中没有切点。';
  }

  /**
   * 统一的图例绘制（自动换行）。
   * 位点线性图谱与质粒环形图谱原本各写了一套图例排版，样式与换行规则还不一致，
   * 这里收敛到一处，两个图谱共用。
   * @param {Array} items [{color, text, shape?}]，shape 为 'circle' 时画圆点，默认方块
   * @returns {number} 图例实际占用的高度
   */
  function drawLegend(ctx, items, x, y, maxW) {
    ctx.font = '10.5px system-ui';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    var cursor = x, row = 0;
    items.forEach(function (it) {
      var w = 16 + ctx.measureText(it.text).width + 14;
      if (cursor + w > x + maxW && cursor > x) { cursor = x; row++; }
      var cy = y + row * 16;
      ctx.fillStyle = it.color;
      if (it.shape === 'circle') {
        ctx.beginPath(); ctx.arc(cursor + 4, cy, 4, 0, Math.PI * 2); ctx.fill();
      } else {
        ctx.fillRect(cursor, cy - 4, 9, 9);
      }
      ctx.fillStyle = COL.label;
      ctx.fillText(it.text, cursor + 13, cy);
      cursor += w;
    });
    return (row + 1) * 16;
  }

  function roundRect(ctx, x, y, w, h, r) {
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.lineTo(x + w - r, y); ctx.quadraticCurveTo(x + w, y, x + w, y + r);
    ctx.lineTo(x + w, y + h - r); ctx.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
    ctx.lineTo(x + r, y + h); ctx.quadraticCurveTo(x, y + h, x, y + h - r);
    ctx.lineTo(x, y + r); ctx.quadraticCurveTo(x, y, x + r, y);
    ctx.closePath();
  }

  function bindMapHover() {
    var cv = $('mapCanvas');
    cv.addEventListener('mousemove', function (ev) {
      if (!mapHits.length) return;
      var rect = cv.getBoundingClientRect();
      var x = ev.clientX - rect.left;
      var best = -1, bestD = 8;
      for (var i = 0; i < mapHits.length; i++) {
        var d0 = Math.abs(mapHits[i].x - x);
        if (d0 < bestD) { bestD = d0; best = i; }
      }
      if (best !== mapHover) { mapHover = best; drawMap(mapData); }
      cv.style.cursor = best >= 0 ? 'pointer' : 'default';
    });
    cv.addEventListener('mouseleave', function () {
      if (mapHover !== -1) { mapHover = -1; drawMap(mapData); }
    });
  }

  /* =========================================================
   * 质粒环形图谱
   * ========================================================= */
  function niceTickStep(n) {
    var raw = n / 10;
    var pow = Math.pow(10, Math.floor(Math.log10(raw)) || 0);
    var cands = [1, 2, 5, 10].map(function (k) { return k * pow; });
    for (var i = 0; i < cands.length; i++) if (cands[i] >= raw) return Math.max(1, Math.round(cands[i]));
    return Math.max(1, Math.round(raw));
  }
  function bpLabel(v) {
    if (v === 0) return '0';
    return v >= 1000 ? num(v / 1000, v % 1000 === 0 ? 0 : 1) + 'k' : String(v);
  }

  function drawPlasmidMap(d) {
    var card = $('plasmidCard');
    var s = state.seq;
    if (!card) return;
    if (!s || !state.circular || !d || state.type === 'protein' || state.type === 'unknown') {
      card.style.display = 'none';
      return;
    }
    card.style.display = '';

    var H = 560;
    var g = setupCanvas($('plasmidCanvas'), H);
    var ctx = g.ctx, W = g.w;
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);

    var cx = W / 2, cy = H / 2 + 4;
    var R = Math.min(W, H) / 2 - 78;
    if (R < 60) R = 60;
    var n = s.length, TAU = Math.PI * 2, START = -Math.PI / 2;  // 从 12 点开始顺时针

    function ang(pos1) { return START + ((pos1 - 1) / n) * TAU; }
    function px(r, a) { return cx + Math.cos(a) * r; }
    function py(r, a) { return cy + Math.sin(a) * r; }

    // 序列圆环
    ctx.strokeStyle = '#cbd5e1'; ctx.lineWidth = 1.6;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, TAU); ctx.stroke();
    ctx.strokeStyle = '#eef2f7'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R - 40, 0, TAU); ctx.stroke();

    // 刻度
    var step = niceTickStep(n);
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.font = '10px system-ui';
    for (var p = 0; p < n; p += step) {
      var a = ang(p + 1);
      ctx.strokeStyle = '#94a3b8'; ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(px(R, a), py(R, a));
      ctx.lineTo(px(R + 6, a), py(R + 6, a));
      ctx.stroke();
      ctx.fillStyle = COL.text;
      ctx.fillText(bpLabel(p), px(R + 18, a), py(R + 18, a));
    }

    // ORF（正链在外圈、负链在内圈）
    var orfs = (n <= 200000) ? getORFs(60).slice(0, 14) : [];
    orfs.forEach(function (o) {
      var r = o.strand === '+' ? R - 22 : R - 33;
      ctx.strokeStyle = o.strand === '+' ? 'rgba(13,148,136,.85)' : 'rgba(124,58,237,.75)';
      ctx.lineWidth = 6;
      ctx.beginPath();
      ctx.arc(cx, cy, r, ang(o.ntStart), ang(o.ntEnd));
      ctx.stroke();
    });

    // 酶切位点
    var enzymes = d.perEnzyme.filter(function (pe) { return pe.hits.length; });
    enzymes.forEach(function (pe, ei) {
      var color = MAP_COLORS[ei % MAP_COLORS.length];
      ctx.strokeStyle = color;
      ctx.lineWidth = 2.2;
      pe.hits.forEach(function (h) {
        var a = ang(h.cut + 1);
        ctx.beginPath();
        ctx.moveTo(px(R - 2, a), py(R - 2, a));
        ctx.lineTo(px(R + 13, a), py(R + 13, a));
        ctx.stroke();
      });
    });

    // 中心信息
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    var name = (state.name || '(未命名)');
    if (name.length > 26) name = name.slice(0, 25) + '…';
    ctx.fillStyle = COL.label; ctx.font = 'bold 13px system-ui';
    ctx.fillText(name, cx, cy - 28);
    ctx.fillStyle = COL.text; ctx.font = '12px system-ui';
    ctx.fillText(n + ' bp · 环状', cx, cy - 8);
    ctx.fillText('GC ' + num(C.gcContent(s), 2) + '%', cx, cy + 10);
    ctx.fillStyle = COL.text; ctx.font = '11px system-ui';
    ctx.fillText('切点 ' + d.perEnzyme.reduce(function (t, pe) { return t + pe.hits.length; }, 0) +
      ' · ORF ' + orfs.length, cx, cy + 28);

    // 图例（底部横向，与其他图谱共用同一套排版）
    var items = enzymes.map(function (pe, ei) {
      return { color: MAP_COLORS[ei % MAP_COLORS.length], text: pe.enzyme.name + ' ×' + pe.hits.length };
    });
    items.push({ color: 'rgba(13,148,136,.85)', text: 'ORF 正链', shape: 'circle' });
    items.push({ color: 'rgba(124,58,237,.75)', text: 'ORF 负链', shape: 'circle' });
    drawLegend(ctx, items, 20, H - 40, W - 40);

    $('plasmidCap').textContent = '环状图谱按序列坐标顺时针排布（12 点方向为第 1 位）。' +
      '径向线段为酶切位点，弧段为 ORF。仅显示最长的 ' + orfs.length + ' 个 ORF（≥60 aa）。';
  }

  /* =========================================================
   * 4. 引物
   * ========================================================= */
  var primerAutoFilled = false;

  function primerOpts() {
    return {
      conc: parseFloat($('prConc').value) || 500,
      na: parseFloat($('prNa').value) || 50,
      mg: parseFloat($('prMg').value) || 0,
      dntp: parseFloat($('prDntp').value) || 0
    };
  }

  function renderPrimer() {
    var fwd = ($('prFwd').value || '').toUpperCase().replace(/\s/g, '');
    var rev = ($('prRev').value || '').toUpperCase().replace(/\s/g, '');
    var opts = primerOpts();

    // 首次打开引物面板时自动从模板两端各取 22 nt 作为示例（用户清空后不再自动填充）
    if (!primerAutoFilled && !fwd && !rev && state.seq && state.seq.length >= 80 &&
        state.type !== 'protein' && state.type !== 'unknown') {
      fwd = state.seq.substr(0, 22);
      rev = C.revComp(state.seq.substr(-22));
      $('prFwd').value = fwd;
      $('prRev').value = rev;
      primerAutoFilled = true;
    }

    if (!fwd && !rev) {
      $('prTbl').querySelector('tbody').innerHTML = '<tr><td colspan="9" style="color:var(--text-3)">输入引物序列，或点「从模板自动挑一对引物」</td></tr>';
      $('prRisk').innerHTML = '';
      $('pcrOut').innerHTML = '<span class="note">—</span>';
      $('pcrProgram').innerHTML = '<span class="note">—</span>';
      $('pcrSeq').textContent = '—';
      return;
    }

    var rows = '';
    [{ n: '正向 (F)', s: fwd }, { n: '反向 (R)', s: rev }].forEach(function (p) {
      if (!p.s) return;
      var tm = C.tmNearestNeighbor(p.s, opts);
      var gc = C.gcContent(p.s);
      var gcT = C.gcClamp(p.s);
      var clampTag = gcT.gc >= 3 ? '<span class="tag ok">良好 ' + gcT.gc + '/5</span>'
        : (gcT.gc >= 2 ? '<span class="tag warn">一般 ' + gcT.gc + '/5</span>' : '<span class="tag bad">偏弱 ' + gcT.gc + '/5</span>');
      var lenTag = (p.s.length >= 18 && p.s.length <= 30) ? '<span class="tag ok">' + p.s.length + ' nt</span>'
        : '<span class="tag warn">' + p.s.length + ' nt</span>';
      var gcTag = (gc >= 40 && gc <= 60) ? '' : ' style="background:#fff7ed"';
      rows += '<tr><td>' + p.n + '</td>' +
        '<td class="num">' + lenTag + '</td>' +
        '<td class="num"' + gcTag + '>' + num(gc, 1) + '</td>' +
        '<td class="num"><b>' + (isNaN(tm) ? '—' : num(tm, 1)) + '</b></td>' +
        '<td class="num">' + num(C.tmGC(p.s), 1) + '</td>' +
        '<td class="num">' + num(C.tmWallace(p.s), 1) + '</td>' +
        '<td class="num">' + num(C.mwNucleic(p.s, 'ssDNA'), 0) + ' Da</td>' +
        '<td style="font-family:var(--mono);font-size:12px">' + esc(gcT.tail) + '</td>' +
        '<td>' + clampTag + '</td></tr>';
    });
    $('prTbl').querySelector('tbody').innerHTML = rows;

    // 风险（基于 ΔG₃₇ 热力学估算，阈值按简化模型标定）
    var risk = '';
    function riskRow(name, val, level, advice) {
      var tag = level === 'ok' ? '<span class="tag ok">低</span>' : (level === 'warn' ? '<span class="tag warn">中</span>' : '<span class="tag bad">高</span>');
      return '<div class="kv"><span class="k">' + name + '</span><span class="v">' + val + ' ' + tag + '</span></div>' +
        (advice ? '<div class="note" style="margin:0 0 4px 12px">' + advice + '</div>' : '');
    }
    function dgText(res) { return 'ΔG = ' + num(res.dG, 2) + ' kcal/mol'; }
    function dgLevel(dG, warnAt, badAt) { return dG <= badAt ? 'bad' : (dG <= warnAt ? 'warn' : 'ok'); }

    if (fwd) {
      var hf = C.hairpinDeltaG(fwd, opts);
      risk += riskRow('正向引物发夹',
        dgText(hf) + (hf.stem ? '（茎 ' + hf.stem + ' bp / 环 ' + hf.loop + ' nt）' : '（未检出茎环）'),
        dgLevel(hf.dG, -1.5, -3.0),
        hf.dG <= -3 ? '发夹结构稳定，会与模板竞争结合导致有效浓度下降，建议调整序列打断茎区。' : '');
      var df = C.dimerDeltaG(fwd, fwd, opts);
      risk += riskRow('正向引物自身二聚体', dgText(df) + (df.run ? '（连续配对 ' + df.run + ' bp）' : ''),
        dgLevel(df.dG, -2.5, -5.0),
        C.dimerTouches3Prime(fwd, fwd, df) && df.dG <= -5 ? '3′ 端参与配对，最易形成引物二聚体并自我扩增，建议改换。' : '');
    }
    if (rev) {
      var hr = C.hairpinDeltaG(rev, opts);
      risk += riskRow('反向引物发夹',
        dgText(hr) + (hr.stem ? '（茎 ' + hr.stem + ' bp / 环 ' + hr.loop + ' nt）' : '（未检出茎环）'),
        dgLevel(hr.dG, -1.5, -3.0), '');
      var dr = C.dimerDeltaG(rev, rev, opts);
      risk += riskRow('反向引物自身二聚体', dgText(dr) + (dr.run ? '（连续配对 ' + dr.run + ' bp）' : ''),
        dgLevel(dr.dG, -2.5, -5.0), '');
    }
    if (fwd && rev) {
      var cross = C.dimerDeltaG(fwd, rev, opts);
      risk += riskRow('交叉二聚体（F × R）', dgText(cross) + (cross.run ? '（连续配对 ' + cross.run + ' bp）' : ''),
        dgLevel(cross.dG, -2.5, -5.0),
        cross.dG <= -5 ? '两条引物互补过强，会互相消耗形成二聚体，建议调整其中一条的 3′ 端。' : '');
      var d3 = C.dimerDeltaG(fwd, rev, opts);
      var touch3 = C.dimerTouches3Prime(fwd, rev, d3);
      var df2 = C.dimerDeltaG(rev, fwd, opts);
      var touch3r = C.dimerTouches3Prime(rev, fwd, df2);
      risk += riskRow('3′ 端交叉配对',
        (touch3 || touch3r) ? '发生在 3′ 端（' + dgText(touch3 ? d3 : df2) + '）' : '未发生在 3′ 端',
        (touch3 || touch3r) ? dgLevel(Math.min(d3.dG, df2.dG), -3.0, -5.0) : 'ok',
        (touch3 || touch3r) ? '3′ 端互补最危险：引物会被聚合酶延伸形成二聚体产物。' : '');

      var tmF = C.tmNearestNeighbor(fwd, opts), tmR = C.tmNearestNeighbor(rev, opts);
      var diff = Math.abs(tmF - tmR);
      risk += riskRow('两引物 Tm 差', num(diff, 1) + ' °C', diff <= 2 ? 'ok' : (diff <= 5 ? 'warn' : 'bad'),
        diff > 3 ? '差值过大时低 Tm 引物会成为限制因素，建议截短高 Tm 引物或延长低 Tm 引物。' : '');
    }
    risk += '<div class="note">ΔG 取 37 °C 自由能，越负越稳定。本模型按连续完美配对累加最近邻自由能并叠加环罚值，' +
      '不处理凸起 / 内部环 / 错配，因此对含错配的稳定结构会低估；阈值已按此标定（发夹 &lt; −3、二聚体 &lt; −5 判为高风险）。</div>';
    $('prRisk').innerHTML = risk;

    // PCR 产物：模板可能有多个引物结合位点，会形成多种产物。
    // 诊断 PCR 里非特异条带常常就是短产物优先扩增 —— 必须全部列出来。
    if (fwd && rev && state.seq && state.type !== 'protein') {
      var prs = C.pcrProducts(state.seq, fwd, rev, { circular: state.circular });
      if (prs.length) {
        var main = prs[0];   // 最短产物最容易被优先扩增
        var outHtml = '';
        if (prs.length === 1) {
          outHtml =
            '<div class="kv"><span class="k">结合位置</span><span class="v">' + main.start + ' – ' + main.end + ' bp</span></div>' +
            '<div class="kv"><span class="k">产物大小</span><span class="v result-ok">' + main.size + ' bp</span></div>' +
            '<div class="kv"><span class="k">占模板比例</span><span class="v">' + num((main.size / state.seq.length) * 100, 1) + '%</span></div>';
        } else {
          outHtml = '<div class="note" style="color:var(--amber)"><b>模板中存在多个引物结合位点</b>，检出 ' +
            prs.length + ' 种可能产物。按长度升序排列：<b>越短的产物越容易被优先扩增</b>，' +
            '若实际电泳出现非特异条带，优先怀疑排在前面的短产物。</div>';
          outHtml += '<div class="tablewrap" style="margin-top:8px"><table class="tbl"><thead><tr>' +
            '<th class="num">#</th><th class="num">产物 (bp)</th><th class="num">区间</th>' +
            '<th class="num">F 位点</th><th class="num">R 位点</th><th></th></tr></thead><tbody>';
          prs.forEach(function (p, i) {
            outHtml += '<tr><td class="num">' + (i + 1) + '</td>' +
              '<td class="num"><b>' + p.size + '</b></td>' +
              '<td class="num">' + p.start + '–' + p.end + '</td>' +
              '<td class="num">#' + p.fwdIndex + '</td>' +
              '<td class="num">#' + p.revIndex + '</td>' +
              '<td><button class="small link" data-pcr="' + i + '">查看</button></td></tr>';
          });
          outHtml += '</tbody></table></div>';
        }
        $('pcrOut').innerHTML = outHtml;
        function showProduct(idx) {
          var p = prs[idx];
          // 跨接缝产物不能简单切片，交给核心层拼接
          var seqTxt = C.pcrProductSeq(state.seq, p);
          $('pcrSeq').textContent = C.formatSeq(seqTxt, 60);
        }
        showProduct(0);
        Array.prototype.forEach.call($('pcrOut').querySelectorAll('button[data-pcr]'), function (b) {
          b.addEventListener('click', function () {
            showProduct(parseInt(b.getAttribute('data-pcr'), 10));
            toast('已显示第 ' + (parseInt(b.getAttribute('data-pcr'), 10) + 1) + ' 种产物序列');
          });
        });
        renderProgram(fwd, rev, main.size, opts);
      } else {
        // pcrProducts 返回空：正向位点、反向位点、或两者相对顺序不成立
        var fHit = state.seq.indexOf(fwd) !== -1;
        var rHit = state.seq.indexOf(C.revComp(rev)) !== -1;
        var why = !fHit ? '模板中未找到正向引物结合位点'
          : (!rHit ? '模板中未找到反向引物结合位点'
            : '反向引物结合位点位于正向引物上游（引物方向可能给反了）');
        $('pcrOut').innerHTML = '<span class="result-bad">未找到产物：' + esc(why) + '</span>';
        $('pcrSeq').textContent = '—';
        $('pcrProgram').innerHTML = '<span class="note">—</span>';
      }
    } else {
      $('pcrOut').innerHTML = '<span class="note">需要模板序列 + 一对引物才能预测产物</span>';
      $('pcrProgram').innerHTML = '<span class="note">—</span>';
    }
  }

  function renderProgram(fwd, rev, size, opts) {
    var tmF = C.tmNearestNeighbor(fwd, opts), tmR = C.tmNearestNeighbor(rev, opts);
    var ta = Math.min(tmF, tmR) - 3;
    var ext = Math.max(15, Math.ceil(size / 1000) * 60);
    var touchdown = ta < 55;
    var rows = [
      ['预变性', '95 °C', '3 min', '1×', '模板充分解链；高 GC 模板可延至 5 min'],
      ['变性', '95 °C', '30 s', '', '也可用 98 °C 10 s（高保真酶）'],
      ['退火', (touchdown ? '65 → ' + ta.toFixed(0) : ta.toFixed(0)) + ' °C', '30 s', '', touchdown ? ' touchdown：每循环降 0.5 °C，共 20 循环后按 ' + ta.toFixed(0) + ' °C 再跑 15 循环' : '取两引物 Tm 较小值 −3 °C'],
      ['延伸', '72 °C', ext + ' s', '', '按 1 kb/min 估算（Taq）；高保真酶约 15–30 s/kb'],
      ['循环数', '', '', '30×', touchdown ? '含 touchdown 20× + 常规 15×' : '质粒模板 25× 通常足够'],
      ['终延伸', '72 °C', '5 min', '1×', '补齐未完成的产物末端'],
      ['保存', '4 °C', '∞', '', '']
    ];
    var html = '<div class="tablewrap"><table class="tbl"><thead><tr><th>步骤</th><th>温度</th><th>时间</th><th class="num">循环</th><th>说明</th></tr></thead><tbody>';
    rows.forEach(function (r) {
      html += '<tr><td><b>' + r[0] + '</b></td><td class="num">' + r[1] + '</td><td class="num">' + r[2] + '</td>' +
        '<td class="num">' + r[3] + '</td><td style="color:var(--text-2)">' + r[4] + '</td></tr>';
    });
    html += '</tbody></table></div>';
    html += '<div class="note">退火温度 ' + ta.toFixed(1) + ' °C 由最近邻 Tm（F ' + num(tmF, 1) + ' / R ' + num(tmR, 1) + ' °C）推算；' +
      '实际最适温度建议用梯度 PCR 在 ' + (ta - 5).toFixed(0) + '–' + (ta + 5).toFixed(0) + ' °C 范围内验证。</div>';
    $('pcrProgram').innerHTML = html;
  }

  /* =========================================================
   * 序列比对
   * ========================================================= */
  var DEMO_ALN_A = 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGTTGCACTTACGCGTACGTACGTAGGCATCCGTAA';
  // 相比 A：中间缺失 GTACGT，并把一处 T 换成 A —— 用于演示 gap 与错配的显示
  var DEMO_ALN_B = 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGTTGCACTTACGCGAAGGCATCCGTAA';

  function doAlign() {
    var a = $('alnA').value, b = $('alnB').value;
    var opts = {
      match: parseInt($('alnMatch').value, 10),
      mismatch: parseInt($('alnMis').value, 10),
      gap: parseInt($('alnGap').value, 10)
    };
    if (isNaN(opts.match)) opts.match = 2;
    if (isNaN(opts.mismatch)) opts.mismatch = -1;
    if (isNaN(opts.gap)) opts.gap = -2;

    // 先判断是否多序列：3 条以上就直接走 MSA，
    // 这样即使 B 框为空（用户把多条序列都贴在 A 框）也能正常比对
    var listA = parseNamedSeqs(a), listB = parseNamedSeqs(b);
    var allSeqs = listA.concat(listB);
    if (allSeqs.length >= 3) {
      renderMSA(allSeqs, opts);
      return;
    }

    if (!a.trim() || !b.trim()) {
      $('alnStats').innerHTML = '<span class="note">请填写两条序列（或把 3 条以上序列贴进 A 框做多序列比对）。</span>';
      $('alnOut').textContent = '—';
      return;
    }

    var res = C.alignPair(a, b, opts);
    if (!res) {
      $('alnStats').innerHTML = '<span class="note">序列为空或只含非法字符。</span>';
      $('alnOut').textContent = '—';
      return;
    }
    if (res.tooBig) {
      $('alnStats').innerHTML = '<span class="result-bad">序列过长：矩阵需要 ' +
        res.cells.toLocaleString() + ' 个单元（上限 ' + res.maxCells.toLocaleString() +
        '）。请把两条序列各截到 2000 nt 以内再比对。</span>';
      $('alnOut').textContent = '—';
      return;
    }

    var mismatches = res.length - res.match - res.gaps;
    var lvl = res.identity >= 99 ? 'result-ok' : (res.identity >= 90 ? '' : 'result-bad');
    $('alnStats').innerHTML =
      '<div class="kv"><span class="k">一致度（identical / 比对长度）</span><span class="v ' + lvl + '">' +
        num(res.identity, 2) + '%</span></div>' +
      '<div class="kv"><span class="k">比对得分</span><span class="v">' + res.score + '</span></div>' +
      '<div class="kv"><span class="k">一致 / 错配 / gap</span><span class="v">' +
        res.match + ' / ' + mismatches + ' / ' + res.gaps + '</span></div>' +
      '<div class="kv"><span class="k">比对长度</span><span class="v">' + res.length + ' (' +
        res.lenA + ' : ' + res.lenB + ')</span></div>';

    // 每 60 列一行，带 1-based 坐标
    var COLS = 60;
    var alnLen = res.alignedA.length;
    var maxRows = 400;
    var lines = [];
    var rows = Math.ceil(alnLen / COLS);
    var limit = Math.min(rows, maxRows);
    for (var r = 0; r < limit; r++) {
      var from = r * COLS, to = Math.min(alnLen, from + COLS);
      var segA = res.alignedA.slice(from, to);
      var segM = res.midline.slice(from, to);
      var segB = res.alignedB.slice(from, to);
      var labelA = pad5(from + 1), labelB = pad5(to);
      lines.push('A ' + labelA + ' ' + segA + ' ' + labelB);
      lines.push('  ' + '      ' + segM);
      lines.push('B ' + labelA + ' ' + segB + ' ' + labelB);
      lines.push('');
    }
    if (rows > maxRows) lines.push('… 共 ' + rows + ' 行，仅显示前 ' + maxRows + ' 行');
    $('alnOut').textContent = lines.join('\n');
  }

  function pad5(v) { var s = String(v); while (s.length < 5) s = ' ' + s; return s; }

  /** 从文本框解析出若干条序列（支持 FASTA 头，也支持一行一条裸序列） */
  function parseNamedSeqs(text) {
    var lines = String(text || '').split(/\r?\n/);
    var out = [], cur = null;
    lines.forEach(function (line) {
      var t = line.trim();
      if (!t) return;
      if (t.charAt(0) === '>') {
        if (cur && cur.seq) out.push(cur);
        cur = { name: t.slice(1).trim(), seq: '' };
      } else {
        if (!cur) cur = { name: '', seq: '' };
        cur.seq += t.replace(/[^A-Za-z]/g, '');
      }
    });
    if (cur && cur.seq) out.push(cur);
    return out.map(function (it) {
      return { name: it.name, seq: it.seq.toUpperCase() };
    }).filter(function (it) { return it.seq.length > 0; });
  }

  /** 多序列比对视图 */
  function renderMSA(list, opts) {
    var msa = C.alignMultiple(list.map(function (it) { return it.seq; }), opts);
    if (!msa) {
      $('alnStats').innerHTML = '<span class="note">至少需要两条有效序列。</span>';
      $('alnOut').textContent = '—';
      return;
    }
    if (msa.tooBig) {
      $('alnStats').innerHTML = '<span class="result-bad">序列过多或过长，无法完成多序列比对' +
        (msa.reason ? '（' + esc(msa.reason) + '）' : '（矩阵 ' + (msa.cells || 0).toLocaleString() + ' 单元）') + '。</span>';
      $('alnOut').textContent = '—';
      return;
    }

    var pairAvg = 0;
    if (msa.pairwise.length) {
      pairAvg = msa.pairwise.reduce(function (a, p) { return a + p.identity; }, 0) / msa.pairwise.length;
    }
    $('alnStats').innerHTML =
      '<div class="kv"><span class="k">序列条数</span><span class="v">' + msa.count + '</span></div>' +
      '<div class="kv"><span class="k">比对列数</span><span class="v">' + msa.columns + '</span></div>' +
      '<div class="kv"><span class="k">完全保守列</span><span class="v result-ok">' + msa.conserved +
        '（' + num(msa.consensusPct, 1) + '%）</span></div>' +
      '<div class="kv"><span class="k">两两一致度均值</span><span class="v">' + num(pairAvg, 2) + '%</span></div>' +
      '<div class="note">中心星算法：以与其他序列平均最相似的一条为参照，其余序列依次并入统一坐标系。' +
      '它速度快、结果稳定，但不等价于 Clustal/MAFFT 那类迭代优化算法，序列长度差异很大时比对质量会下降。</div>';

    // 保守列标记行
    var consLine = '';
    for (var c = 0; c < msa.columns; c++) {
      var first = msa.rows[0].charAt(c), same = true;
      for (var k = 1; k < msa.rows.length; k++) {
        if (msa.rows[k].charAt(c) !== first) { same = false; break; }
      }
      consLine += (same && first !== '-') ? '|' : ' ';
    }

    var COLS = 60, lines = [], width = msa.columns;
    var labelW = 0;
    list.forEach(function (it, i) {
      var n = (it.name || ('#' + (i + 1)));
      labelW = Math.max(labelW, Math.min(n.length, 18));
    });
    function lab(s) { var t = String(s); if (t.length > 18) t = t.slice(0, 17) + '…'; while (t.length < labelW) t += ' '; return t; }

    for (var from = 0; from < width; from += COLS) {
      var to = Math.min(width, from + COLS);
      msa.rows.forEach(function (row, i) {
        lines.push(lab(list[i] ? (list[i].name || ('#' + (i + 1))) : ('#' + (i + 1))) + ' ' + row.slice(from, to));
        if (i === 0) lines.push(new Array(labelW + 2).join(' ') + consLine.slice(from, to));
      });
      lines.push('');
    }
    $('alnOut').textContent = lines.join('\n');
  }

  /* =========================================================
   * 虚拟克隆
   * ========================================================= */
  function cleanSeqText(t) {
    return String(t || '').toUpperCase().replace(/[^ACGTUNRYKMSWBDHV]/g, '');
  }

  /** 把 51 种酶填进下拉框 */
  function fillEnzymeSelect(sel, allowNone) {
    if (!sel) return;
    var cur = sel.value;
    sel.innerHTML = '';
    if (allowNone) {
      var none = document.createElement('option');
      none.value = '';
      none.textContent = '（不酶切，用完整序列）';
      sel.appendChild(none);
    }
    E.list.forEach(function (en) {
      var o = document.createElement('option');
      o.value = en.name;
      o.textContent = en.name + '  ' + en.site;
      sel.appendChild(o);
    });
    if (cur) sel.value = cur;
  }

  function cloneFragments(side) {
    var isBack = side === 'back';
    var seq = cleanSeqText($(isBack ? 'clBackSeq' : 'clInsSeq').value);
    var enzName = $(isBack ? 'clBackEnz' : 'clInsEnz').value;
    if (!seq) return { seq: '', frags: [] };
    if (!enzName) {
      // 不酶切：整条作为一个片段
      return {
        seq: seq,
        frags: [{ seq: seq, length: seq.length, leftEnd: null, rightEnd: null, uncut: true, index: 0 }]
      };
    }
    var en = E.byName(enzName);
    if (!en) return { seq: seq, frags: [] };
    // 载体按环状、插入按线性
    var frags = C.digestOneDetailed(seq, en, !!isBack);
    frags.forEach(function (f, i) { f.index = i; });
    return { seq: seq, frags: frags };
  }

  function endText(e) {
    if (!e) return '<span class="tag">线性末端（无需配对）</span>';
    if (e.type === 'blunt') return '<span class="tag ok">平端</span>';
    if (e.degenerate) return '<span class="tag warn">简并突出端（无法判定）</span>';
    return '<span class="tag ' + (e.type === '5' ? 'brand' : '') + '">' + e.type + '′ 突出 ' + e.overhang + '</span>';
  }

  /** 刷新某一侧的片段下拉框 */
  function refreshCloneFrags(side) {
    var isBack = side === 'back';
    var sel = $(isBack ? 'clBackFrag' : 'clInsFrag');
    var info = cloneFragments(side);
    var cur = sel.value;
    sel.innerHTML = '';
    if (!info.seq) {
      var o0 = document.createElement('option');
      o0.value = ''; o0.textContent = '（请先填入序列）';
      sel.appendChild(o0);
      return;
    }
    if (!info.frags.length) {
      var o1 = document.createElement('option');
      o1.value = ''; o1.textContent = '（没有可用片段）';
      sel.appendChild(o1);
      return;
    }
    info.frags.forEach(function (f, i) {
      var o = document.createElement('option');
      o.value = String(i);
      o.textContent = '#' + (i + 1) + '  ' + f.length + ' bp' +
        (f.wrap ? '（跨接缝）' : '') + (f.uncut ? '（未切开）' : '');
      sel.appendChild(o);
    });
    // 默认选最长的那段（通常是骨架/目标片段）
    var longest = 0;
    info.frags.forEach(function (f, i) { if (f.length > info.frags[longest].length) longest = i; });
    sel.value = (info.frags.length > 1) ? String(longest) : '0';
    if (cur && info.frags[parseInt(cur, 10)]) sel.value = cur;
  }

  function doClone() {
    var bInfo = cloneFragments('back');
    var iInfo = cloneFragments('ins');
    if (!bInfo.seq || !iInfo.seq) {
      $('clOut').innerHTML = '<span class="note">请先填入载体与插入片段的序列。</span>';
      $('clSeq').textContent = '—';
      return;
    }
    var bf = bInfo.frags[parseInt($('clBackFrag').value, 10)];
    var inf = iInfo.frags[parseInt($('clInsFrag').value, 10)];
    if (!bf || !inf) {
      $('clOut').innerHTML = '<span class="note">请选择载体骨架片段与插入片段。</span>';
      $('clSeq').textContent = '—';
      return;
    }

    var res = C.ligateInsert(
      { seq: bf.seq, leftEnd: bf.leftEnd, rightEnd: bf.rightEnd },
      { seq: inf.seq, leftEnd: inf.leftEnd, rightEnd: inf.rightEnd }
    );

    var out = '';
    out += '<div class="kv"><span class="k">载体骨架</span><span class="v">' + bf.length + ' bp' +
      (bf.uncut ? '（未酶切）' : '（切点 ' + bf.leftCut + ' / ' + bf.rightCut + '）') + '</span></div>';
    out += '<div class="kv"><span class="k">插入片段</span><span class="v">' + inf.length + ' bp' +
      (inf.uncut ? '（未酶切）' : '（切点 ' + inf.leftCut + ' / ' + inf.rightCut + '）') + '</span></div>';

    res.junctions.forEach(function (j) {
      var tag = j.ok ? '<span class="tag ok">可连接</span>' : '<span class="tag bad">不可连接</span>';
      out += '<div class="kv"><span class="k">' + j.name + '</span><span class="v">' +
        endText(j.byBackbone) + ' ⟷ ' + endText(j.byInsert) + ' ' + tag + '</span></div>';
    });

    if (res.ok) {
      out += '<div class="kv"><span class="k">重组质粒</span><span class="v result-ok">' +
        res.size + ' bp（相对载体 ' + (res.size - bInfo.seq.length >= 0 ? '+' : '') +
        (res.size - bInfo.seq.length) + ' bp）</span></div>';
      out += '<div class="kv"><span class="k">GC 含量</span><span class="v">' + num(C.gcContent(res.seq), 2) + '%</span></div>';
    } else {
      out += '<div class="note" style="color:var(--amber)">' + esc(res.reason) + '</div>';
      out += '<div class="note">提示：同尾酶（如 BamHI/BglII、SalI/XhoI）产生的黏端可以互相连接；' +
        '若两端不兼容，可考虑补平后平端连接，或改用带相同酶切位点的引物重新扩增。</div>';
    }
    $('clOut').innerHTML = out;
    $('clSeq').textContent = C.formatSeq(res.seq, 60);
    window.__cloneSeq = res.seq;
  }

  /* =========================================================
   * 4.5 蛋白 & 密码子
   * ========================================================= */
  var AA_CLASS = {
    A: '非极性', V: '非极性', I: '非极性', L: '非极性', M: '非极性', P: '非极性',
    G: '极性', S: '极性', T: '极性', N: '极性', Q: '极性', C: '极性（含硫）',
    F: '芳香族', W: '芳香族', Y: '芳香族',
    D: '酸性', E: '酸性', K: '碱性', R: '碱性', H: '碱性'
  };

  function readProtSeq() {
    return ($('protSeq').value || '')
      .replace(/^>[^\n]*(\n|$)/gm, '')
      .replace(/[^A-Za-z*]/g, '').toUpperCase();
  }

  function analyzeProtein() {
    var seq = readProtSeq().replace(/\*/g, '');
    var pp = C.proteinProperties(seq);
    if (!pp) {
      $('protStats').innerHTML = '<span class="note">请先提供蛋白序列（可点上方「提取最长 ORF」）。</span>';
      $('protCompTbl').querySelector('tbody').innerHTML = '';
      $('protNote').textContent = '';
      return;
    }
    var tiles = [];
    function tile(k, v, u) { tiles.push('<div class="stat"><div class="k">' + k + '</div><div class="v">' + v + (u ? '<span class="u">' + u + '</span>' : '') + '</div></div>'); }
    tile('长度', pp.length, 'aa');
    tile('分子量', pp.mw >= 10000 ? num(pp.mw / 1000, 2) : num(pp.mw, 1), pp.mw >= 10000 ? 'kDa' : 'Da');
    tile('等电点 pI', num(pp.pI, 2));
    tile('pH 7 净电荷', (pp.charge7 >= 0 ? '+' : '') + num(pp.charge7, 1));
    tile('GRAVY 疏水性', num(pp.gravy, 3));
    tile('脂肪族指数', num(pp.aliphatic, 1));
    tile('不稳定指数', num(pp.instability, 1), pp.unstable ? '⚠' : '✓');
    tile('ε₂₈₀（全氧化）', num(pp.extinctionOx / 1000, 1), 'mM⁻¹cm⁻¹');
    tile('A₂₈₀ 0.1%（氧化）', num(pp.abs01Ox, 3));
    $('protStats').innerHTML = tiles.join('');

    var tb = $('protCompTbl').querySelector('tbody');
    var rows = [];
    'ACDEFGHIKLMNPQRSTVWY'.split('').forEach(function (r) {
      var n = pp.counts[r] || 0;
      rows.push('<tr><td><b>' + r + '</b></td><td class="num">' + n + '</td><td class="num">' +
        num(100 * n / pp.length, 1) + '%</td><td>' + AA_CLASS[r] + '</td></tr>');
    });
    tb.innerHTML = rows.join('');

    var notes = [];
    notes.push(pp.unstable
      ? '不稳定指数 ' + num(pp.instability, 1) + ' > 40：该蛋白在体外可能不稳定，纯化时建议加保护剂并低温操作。'
      : '不稳定指数 ' + num(pp.instability, 1) + ' ≤ 40：预测为稳定蛋白。');
    if (pp.gravy > 0.5) notes.push('GRAVY 明显为正：疏水性强，可能是膜蛋白或需要促溶标签。');
    if (pp.extinctionOx === 0) notes.push('不含 W/Y：无法用 A₂₈₀ 定量，需改用 BCA/Bradford。');
    $('protNote').textContent = notes.join(' ');
  }

  function codonOptimize() {
    var seq = readProtSeq();
    if (!seq.replace(/\*/g, '')) {
      $('codonOut').innerHTML = '<span class="note">请先在上方提供蛋白序列。</span>';
      $('codonDna').textContent = '—';
      return;
    }
    var rt = C.reverseTranslate(seq);
    var out = '';
    out += '<div class="kv"><span class="k">DNA 长度</span><span class="v">' + rt.dna.length + ' bp（' +
      seq.replace(/\*/g, '').length + ' aa' + (rt.skipped ? '，跳过非法残基 ' + rt.skipped + ' 个' : '') + '）</span></div>';
    out += '<div class="kv"><span class="k">GC 含量</span><span class="v">' + num(C.gcContent(rt.dna), 1) + '%</span></div>';
    out += '<div class="kv"><span class="k">回译校验</span><span class="v">' +
      (rt.backTranslateOk ? '<span class="tag ok">与原蛋白完全一致</span>' : '<span class="tag bad">不一致，请检查输入</span>') + '</span></div>';
    // 新引入的常用酶切位点（克隆载体设计时通常要避开）
    var hits = [];
    E.list.forEach(function (en) {
      var h = C.findSites(rt.dna, en, false);
      if (h.length) hits.push({ name: en.name, n: h.length });
    });
    if (hits.length) {
      out += '<div class="kv"><span class="k">产物内酶切位点</span><span class="v">' +
        hits.map(function (h) { return '<span class="tag warn">' + esc(h.name) + ' ×' + h.n + '</span>'; }).join(' ') + '</span></div>';
      out += '<div class="note">提示：若计划用上述酶做亚克隆，可在对应位置改用同义密码子消灭位点（当前版本按最高频密码子统一选码，尚未逐位点规避）。</div>';
    } else {
      out += '<div class="kv"><span class="k">产物内酶切位点</span><span class="v"><span class="tag ok">51 种常用酶均无切点</span></span></div>';
    }
    $('codonOut').innerHTML = out;
    $('codonDna').textContent = C.formatSeq(rt.dna, 60);
    window.__codonDna = rt.dna;
  }

  function bindProtein() {
    $('btnProtFromOrf').addEventListener('click', function () {
      if (!state.seq || state.type === 'protein') { toast('当前序列不是核酸，无法提取 ORF'); return; }
      var orfs = getORFs(30);
      if (!orfs.length) { toast('未找到 ≥30 aa 的 ORF'); return; }
      $('protSeq').value = '>longest_orf ' + orfs[0].aaLength + ' aa\n' + orfs[0].protein;
      toast('已提取最长 ORF（' + orfs[0].aaLength + ' aa）');
      analyzeProtein();
    });
    $('btnProtFromSeq').addEventListener('click', function () {
      if (!state.seq) { toast('没有序列'); return; }
      $('protSeq').value = '>' + (state.name || 'protein') + '\n' + state.seq;
      analyzeProtein();
    });
    $('btnProtClear').addEventListener('click', function () { $('protSeq').value = ''; });
    $('btnProtAnalyze').addEventListener('click', analyzeProtein);
    $('btnCodonOpt').addEventListener('click', codonOptimize);
    $('btnCodonCopy').addEventListener('click', function () {
      if (!window.__codonDna) { toast('请先生成'); return; }
      copyText(window.__codonDna, '优化 DNA');
    });
    $('btnCodonToInput').addEventListener('click', function () {
      if (!window.__codonDna) { toast('请先生成'); return; }
      $('seqInput').value = '>codon_optimized_ecoli\n' + C.formatSeq(window.__codonDna, 60);
      $('chkCircular').checked = false;
      readInput();
      toast('已载入为当前序列');
    });
  }

  /* =========================================================
   * 4.6 Motif 搜索 & 双酶切速查
   * ========================================================= */
  function doMotif() {
    var tb = $('motifTbl').querySelector('tbody');
    if (!state.seq || state.type === 'protein') {
      tb.innerHTML = ''; $('motifOut').textContent = '请先在上方输入核酸序列。';
      return;
    }
    var pat = ($('motifPat').value || '').toUpperCase().replace(/\s/g, '');
    if (!pat) { tb.innerHTML = ''; $('motifOut').textContent = '请输入模式（支持 IUPAC 简并碱基）。'; return; }
    if (!/^[ACGTURYSWKMBDHVN]+$/.test(pat)) {
      tb.innerHTML = ''; $('motifOut').textContent = '模式含无法识别的字符：仅允许 A C G T U R Y S W K M B D H V N。';
      return;
    }
    pat = pat.replace(/U/g, 'T');
    var MAX = 500;
    var hits = C.findMotifs(state.seq, pat, {
      circular: state.circular,
      bothStrands: $('motifBoth').checked
    });
    var rows = hits.slice(0, MAX).map(function (h, i) {
      return '<tr><td class="num">' + (i + 1) + '</td><td class="num">' + h.pos +
        (h.wrap ? ' <span class="tag warn">跨接缝</span>' : '') + '</td><td>' +
        (h.strand === '+' ? '<span class="tag brand">正链</span>' : '<span class="tag">负链</span>') +
        '</td><td class="num" style="text-align:left;font-family:var(--mono)">' + esc(h.match) + '</td></tr>';
    });
    tb.innerHTML = rows.join('') || '<tr><td colspan="4" style="text-align:center;color:var(--text-3)">无命中</td></tr>';
    $('motifOut').textContent = '共 ' + hits.length + ' 处命中' +
      (hits.truncated ? '（已达上限并截断，请缩短序列或提高模式特异性）' : '') +
      (hits.length > MAX ? '，仅显示前 ' + MAX + ' 条' : '') +
      '（' + (state.circular ? '环状' : '线性') + ($('motifBoth').checked ? '，双链' : '，仅正链') + '）。';
  }

  function doDD() {
    if (!state.seq || state.type === 'protein') {
      $('ddOut').innerHTML = '<span class="note">请先在上方输入核酸序列。</span>';
      return;
    }
    var a = E.byName($('ddEnzA').value), b = E.byName($('ddEnzB').value);
    if (!a || !b) { $('ddOut').innerHTML = '<span class="note">请选择两种酶。</span>'; return; }
    var r = C.doubleDigestSummary(state.seq, a, b, state.circular);
    var out = '';
    out += '<div class="kv"><span class="k">' + esc(a.name) + '</span><span class="v">' +
      r.cutsA + ' 个切点 · ' + endText(r.endsA) + '</span></div>';
    out += '<div class="kv"><span class="k">' + esc(b.name) + '</span><span class="v">' +
      r.cutsB + ' 个切点 · ' + endText(r.endsB) + '</span></div>';
    out += '<div class="kv"><span class="k">黏端互补（同尾）</span><span class="v">' +
      (r.compatible ? '<span class="tag warn">互补，可互连</span>' : '<span class="tag ok">不互补</span>') + '</span></div>';
    out += '<div class="kv"><span class="k">联合酶切片段</span><span class="v">' +
      (r.fragments.length ? r.fragments.join(' / ') + ' bp' : '（无切点）') + '</span></div>';
    r.notes.forEach(function (n) { out += '<div class="note">· ' + esc(n) + '</div>'; });
    $('ddOut').innerHTML = out;
  }

  function bindDigestExtra() {
    fillEnzymeSelect($('ddEnzA'));
    fillEnzymeSelect($('ddEnzB'));
    $('ddEnzA').value = 'EcoRI';
    $('ddEnzB').value = 'HindIII';
    $('btnDD').addEventListener('click', doDD);
    $('motifPreset').addEventListener('change', function () {
      if ($('motifPreset').value) $('motifPat').value = $('motifPreset').value;
      $('motifPreset').value = '';
    });
    $('btnMotif').addEventListener('click', doMotif);
  }

  /* =========================================================
   * 4.7 定点突变引物
   * ========================================================= */
  function doMut() {
    if (!state.seq || state.type === 'protein') {
      $('mutOut').innerHTML = '<span class="note">请先在上方输入模板核酸序列（质粒请勾选「环状」）。</span>';
      window.__mutRes = null;
      return;
    }
    var pos = parseInt($('mutPos').value, 10);
    var res = C.designMutPrimers(state.seq, pos, $('mutNew').value, {
      oldLen: parseInt($('mutOld').value, 10) || 0,
      targetTm: parseFloat($('mutTm').value) || 62,
      circular: state.circular,
      dnaConc: parseFloat($('prConc').value) || 500,
      naConc: parseFloat($('prNa').value) || 50
    });
    if (!res.ok) {
      $('mutOut').innerHTML = '<span class="note" style="color:var(--red)">' + esc(res.reason) + '</span>';
      window.__mutRes = null;
      return;
    }
    window.__mutRes = res;
    var out = '';
    out += '<div class="kv"><span class="k">突变</span><span class="v">' + esc(res.mutDesc) + '</span></div>';
    out += '<div class="kv"><span class="k">正向引物（5′→3′）</span><span class="v" style="font-family:var(--mono)">' +
      '<b style="color:var(--red)">' + esc(res.fwdTail || '') + '</b>' + esc(res.fwdAnneal) +
      ' <span class="tag">' + res.fwd.length + ' nt</span></span></div>';
    out += '<div class="kv"><span class="k">反向引物（5′→3′）</span><span class="v" style="font-family:var(--mono)">' +
      esc(res.rev) + ' <span class="tag">' + res.rev.length + ' nt</span></span></div>';
    out += '<div class="kv"><span class="k">退火区 Tm（最近邻）</span><span class="v">F ' + num(res.tmF, 1) +
      ' °C / R ' + num(res.tmR, 1) + ' °C（红字为突变碱基，不计 Tm）</span></div>';
    out += '<div class="kv"><span class="k">产物（突变质粒）</span><span class="v">' + res.productLen + ' bp</span></div>';
    out += '<div class="note">程序建议（Q5 类高保真酶）：98 °C 30 s；98 °C 10 s → ' +
      num(Math.min(res.tmF, res.tmR), 0) + ' °C 20 s → 72 °C ' +
      Math.max(15, Math.round(res.productLen / 1000 * 30)) + ' s，25 循环；72 °C 2 min。产物经 KLD 处理（激酶-连接酶-DpnI）后转化。</div>';
    res.warnings.forEach(function (w) {
      out += '<div class="note" style="color:var(--amber)">⚠ ' + esc(w) + '</div>';
    });
    $('mutOut').innerHTML = out;
  }

  function bindMut() {
    $('btnMut').addEventListener('click', doMut);
    $('btnMutToPrimers').addEventListener('click', function () {
      var r = window.__mutRes;
      if (!r || !r.ok) { toast('请先成功设计突变引物'); return; }
      $('prFwd').value = r.fwd;
      $('prRev').value = r.rev;
      toast('已填入引物输入框，可点「分析引物」做完整评估');
    });
  }

  /* =========================================================
   * 4.8 连接反应用量
   * ========================================================= */
  function doLig() {
    var r = C.ligationAmounts($('ligVecBp').value, $('ligInsBp').value, $('ligVecNg').value, $('ligRatio').value);
    if (!r) {
      $('ligOut').innerHTML = '<span class="note" style="color:var(--red)">输入不完整：载体/插入长度、载体用量与摩尔比都必须大于 0。</span>';
      return;
    }
    var out = '';
    out += '<div class="kv"><span class="k">载体</span><span class="v">' + $('ligVecNg').value + ' ng = ' +
      num(r.vectorPmol * 1000, 1) + ' fmol</span></div>';
    out += '<div class="kv"><span class="k">插入片段（' + $('ligRatio').value + ' : 1）</span><span class="v result-ok">' +
      num(r.insertNg, 1) + ' ng（' + num(r.insertPmol * 1000, 1) + ' fmol）</span></div>';
    var vc = parseFloat($('ligVecConc').value), ic = parseFloat($('ligInsConc').value);
    if (vc > 0) out += '<div class="kv"><span class="k">载体吸取体积</span><span class="v">' + num(parseFloat($('ligVecNg').value) / vc, 2) + ' µL</span></div>';
    if (ic > 0) out += '<div class="kv"><span class="k">插入吸取体积</span><span class="v">' + num(r.insertNg / ic, 2) + ' µL</span></div>';
    if (!(vc > 0) || !(ic > 0)) out += '<div class="note">填入两侧浓度后可显示各需吸取的 µL 数。</div>';
    out += '<div class="note">常规 T4 连接（20 µL 体系）：载体 + 插入 + 2 µL 10× 连接缓冲液 + 1 µL T4 连接酶，补水至 20 µL，16 °C 过夜或室温 10–30 min（黏端）。</div>';
    $('ligOut').innerHTML = out;
  }

  function bindLig() {
    $('btnLig').addEventListener('click', doLig);
    $('btnLigFill').addEventListener('click', function () {
      var bInfo = cloneFragments('back');
      var iInfo = cloneFragments('ins');
      var bf = bInfo.frags[parseInt($('clBackFrag').value, 10)];
      var inf = iInfo.frags[parseInt($('clInsFrag').value, 10)];
      if (!bf || !inf) { toast('请先在上方选好载体骨架与插入片段'); return; }
      $('ligVecBp').value = bf.length;
      $('ligInsBp').value = inf.length;
      doLig();
      toast('已按所选片段填入长度（骨架 ' + bf.length + ' bp / 插入 ' + inf.length + ' bp）');
    });
  }

  /* =========================================================
   * 5. 计算器
   * ========================================================= */
  function bindCalc() {
    $('btnDilute').addEventListener('click', function () {
      var f = { C1: $('dC1'), V1: $('dV1'), C2: $('dC2'), V2: $('dV2') };
      var vals = {}, empties = [];
      Object.keys(f).forEach(function (k) {
        var v = f[k].value.trim();
        vals[k] = v === '' ? null : parseFloat(v);
        if (vals[k] === null) empties.push(k);
      });
      if (empties.length !== 1) {
        $('diluteOut').innerHTML = '<b style="color:var(--amber)">请恰好留空一项</b>（当前留空：' + (empties.join(', ') || '无') + '）';
        return;
      }
      var k = empties[0], res;
      if (k === 'V1') res = (vals.C2 * vals.V2) / vals.C1;
      else if (k === 'C1') res = (vals.C2 * vals.V2) / vals.V1;
      else if (k === 'V2') res = (vals.C1 * vals.V1) / vals.C2;
      else res = (vals.C1 * vals.V1) / vals.V2;
      f[k].value = Number(res.toFixed(4));
      var unit = { C1: '', V1: '', C2: '', V2: '' };
      $('diluteOut').innerHTML = '<b class="result-ok">' + k + ' = ' + Number(res.toFixed(4)) + '</b>' +
        (k === 'V1' ? '（取母液体积；若需配制，另加 ' + Number((vals.V2 - res).toFixed(4)) + ' 体积单位的稀释液）' : '');
    });

    $('btnMolFromMass').addEventListener('click', function () {
      var mw = parseFloat($('mMW').value), m = parseFloat($('mMass').value);
      if (!mw || isNaN(m)) { $('molOut').textContent = '请填写分子量与质量'; return; }
      $('mMol').value = Number(((m / mw) * 1000).toFixed(6));
      $('molOut').innerHTML = m + ' mg ÷ ' + mw + ' g/mol = <b class="result-ok">' + ((m / mw) * 1000).toFixed(4) + ' mmol</b> = ' + ((m / mw) * 1e6).toFixed(1) + ' nmol';
    });
    $('btnMassFromMol').addEventListener('click', function () {
      var mw = parseFloat($('mMW').value), n = parseFloat($('mMol').value);
      if (!mw || isNaN(n)) { $('molOut').textContent = '请填写分子量与物质的量'; return; }
      $('mMass').value = Number((n * mw).toFixed(4));
      $('molOut').innerHTML = n + ' mmol × ' + mw + ' g/mol = <b class="result-ok">' + (n * mw).toFixed(4) + ' mg</b>';
    });
    $('btnPct').addEventListener('click', function () {
      var pctv = parseFloat($('pPct').value), vol = parseFloat($('pVol').value), mw = parseFloat($('pMW').value);
      if (isNaN(pctv) || isNaN(vol) || !mw) { $('molOut').textContent = '请填写完整'; return; }
      var massG = (pctv / 100) * vol;                // g
      var molar = (massG / mw) / (vol / 1000);        // mol/L
      $('molOut').innerHTML = pctv + '% (w/v) × ' + vol + ' mL → 称取 <b class="result-ok">' + (massG * 1000).toFixed(2) + ' mg</b>；' +
        '摩尔浓度 = <b>' + (molar * 1000).toFixed(1) + ' mM</b> = ' + (molar * 1e6).toFixed(0) + ' µM';
    });

    $('btnRx').addEventListener('click', function () {
      var V = parseFloat($('rxVol').value), N = parseFloat($('rxN').value), EXn = parseFloat($('rxExtra').value) || 0;
      var tplC = parseFloat($('rxTpl').value), tplNg = parseFloat($('rxTplNg').value);
      var pStock = parseFloat($('rxPrStock').value), pFinal = parseFloat($('rxPrFinal').value);
      var dStock = parseFloat($('rxDntpStock').value), dFinal = parseFloat($('rxDntpFinal').value);
      var bStock = parseFloat($('rxBufStock').value);
      var eStock = parseFloat($('rxPolStock').value), eU = parseFloat($('rxPolU').value);
      if ([V, N, tplC, pStock, dStock, bStock, eStock].some(isNaN)) { $('rxOut').textContent = '参数填写不完整'; return; }

      var nMix = N + EXn;
      var items = [];
      items.push({ name: '模板 DNA', per: tplNg / tplC, mix: false, final: tplNg + ' ng/管', note: tplC + ' ng/µL 母液' });
      var pVol = (pFinal / 1000) / pStock * V;   // nM→µM
      items.push({ name: '正向引物', per: pVol, final: pFinal + ' nM', note: pStock + ' µM 母液' });
      items.push({ name: '反向引物', per: pVol, final: pFinal + ' nM', note: pStock + ' µM 母液' });
      var dVol = (dFinal / 1000) / dStock * V;
      items.push({ name: 'dNTP Mix', per: dVol, final: dFinal + ' µM each', note: dStock + ' mM each 母液' });
      items.push({ name: '缓冲液 ' + bStock + '×', per: V / bStock, final: '1×', note: '' });
      var eVol = eU / eStock;
      items.push({ name: 'DNA 聚合酶', per: eVol, final: eU + ' U/管', note: eStock + ' U/µL' });

      var Water = V - items.reduce(function (a, b) { return a + b.per; }, 0);
      if (Water < 0) {
        items.push({ name: 'ddH₂O', per: 0, final: '', note: '体积已超总体系' });
      } else {
        items.push({ name: 'ddH₂O', per: Water, mix: true, final: '', note: '补足至 ' + V + ' µL' });
      }

      var html = '';
      var sumPer = 0, sumMix = 0;
      items.forEach(function (it) {
        var isTpl = it.name === '模板 DNA';
        // 模板与酶通常每管单独加，不计入 Master Mix
        var mixVol = (isTpl || it.name === 'DNA 聚合酶') ? 0 : it.per * nMix;
        sumPer += it.per;
        sumMix += mixVol;
        html += '<tr><td>' + it.name + '</td>' +
          '<td class="num">' + num(it.per, 2) + '</td>' +
          '<td class="num">' + (mixVol > 0 ? num(mixVol, 2) : '<span style="color:var(--text-3)">单独加</span>') + '</td>' +
          '<td class="num">' + (it.final || '—') + '</td>' +
          '<td style="color:var(--text-2)">' + esc(it.note) + '</td></tr>';
      });
      html += '</tbody>';
      $('rxTbl').querySelector('tbody').innerHTML = html;
      $('rxTbl').querySelector('tfoot').innerHTML =
        '<tr><th>合计</th><th class="num">' + num(sumPer, 2) + ' µL</th>' +
        '<th class="num">' + num(sumMix, 2) + ' µL</th><th colspan="2"></th></tr>';

      var warn = Water < 0
        ? '<b class="result-bad">组分体积之和已超过 ' + V + ' µL，请降低模板量或提高母液浓度。</b>'
        : '按 <b>' + nMix + '</b> 管（' + N + ' 管 + ' + EXn + ' 管余量）配制 Master Mix；每管分装 ' +
          num(sumMix / nMix, 2) + ' µL 后，再各加模板 ' + num(tplNg / tplC, 2) + ' µL 与酶 ' + num(eU / eStock, 2) + ' µL。';
      $('rxOut').innerHTML = warn;
    });

    function rpm2g() {
      var r = parseFloat($('cfR').value), rpm = parseFloat($('cfRpm').value);
      if (!r || isNaN(rpm)) { $('cfOut').textContent = '请填写半径与转速'; return; }
      var g = 1.118e-5 * r * rpm * rpm;
      $('cfG').value = Number(g.toFixed(0));
      $('cfOut').innerHTML = r + ' cm · ' + rpm + ' rpm → <b class="result-ok">' + g.toFixed(0) + ' ×g</b>';
    }
    function g2rpm() {
      var r = parseFloat($('cfR').value), g = parseFloat($('cfG').value);
      if (!r || isNaN(g)) { $('cfOut').textContent = '请填写半径与离心力'; return; }
      var rpm = Math.sqrt(g / (1.118e-5 * r));
      $('cfRpm').value = Number(rpm.toFixed(0));
      $('cfOut').innerHTML = g + ' ×g @ ' + r + ' cm → <b class="result-ok">' + rpm.toFixed(0) + ' rpm</b>';
    }
    $('btnRpm2g').addEventListener('click', rpm2g);
    $('btnG2rpm').addEventListener('click', g2rpm);

    $('btnOd6').addEventListener('click', function () {
      var od = parseFloat($('od6').value), dil = parseFloat($('od6dil').value) || 1;
      var k = parseFloat($('od6kind').value);
      if (isNaN(od)) { $('od6Out').textContent = '请填写 OD₆₀₀'; return; }
      var cells = od * dil * k;
      $('od6Out').innerHTML = 'OD₆₀₀ = ' + od + ' × ' + dil + ' 稀释 → 约 <b class="result-ok">' +
        cells.toExponential(2) + ' cells/mL</b>（' + (cells / 1000).toExponential(2) + ' cells/µL）。' +
        '经验系数仅供估算，精确计数请用血球计数板或流式。';
    });
  }

  /* =========================================================
   * 事件绑定
   * ========================================================= */
  function bind() {
    // 首次进入时预选一对经典克隆用酶（双酶切），避免酶切面板空空如也或片段过碎
    state.selectedEnzymes = ['EcoRI', 'HindIII'];

    var inputTimer = null;
    $('seqInput').addEventListener('input', function () {
      clearTimeout(inputTimer);
      // 自适应防抖：大序列一次完整渲染要几百毫秒，编辑过程中没必要每停顿一次就跑。
      // 小序列保持 220 ms 的即时反馈，超过 50 kb 放宽到 500 ms。
      var len = $('seqInput').value.length;
      inputTimer = setTimeout(readInput, len > 50000 ? 500 : 220);
    });
    $('chkCircular').addEventListener('change', readInput);

    // 示例
    var sel = $('exampleSel');
    EX.forEach(function (e, i) {
      var o = document.createElement('option');
      o.value = String(i);
      o.textContent = e.name + '（' + e.seq.length + ' bp）';
      sel.appendChild(o);
    });
    sel.addEventListener('change', function () {
      if (sel.value === '') return;
      var e = EX[parseInt(sel.value, 10)];
      $('seqInput').value = '>' + e.id + ' ' + e.name + '\n' + C.formatSeq(e.seq, 60);
      $('chkCircular').checked = !!e.circular;
      readInput();
      toast('已载入 ' + e.name);
    });

    $('histSel').addEventListener('change', function () {
      var sel = $('histSel');
      if (sel.value === '') return;
      var it = histLoad()[parseInt(sel.value, 10)];
      if (!it || !it.seq) { toast('该条历史已失效'); renderHistory(); return; }
      $('seqInput').value = '>' + (it.name || 'history') + '\n' + C.formatSeq(it.seq, 60);
      $('chkCircular').checked = !!it.circular;
      readInput();
      sel.value = '';
      toast('已载入历史序列（' + it.seq.length + ' nt）');
    });
    $('btnHistClear').addEventListener('click', function () {
      histClear();
      renderHistory();
      toast('已清空本机历史记录');
    });

    $('btnDemoRandom').addEventListener('click', function () {
      var s = C.randomSeq(2000, 50);
      $('seqInput').value = '>random_2000bp\n' + C.formatSeq(s, 60);
      $('chkCircular').checked = false;
      readInput();
      toast('已生成 2000 bp 随机序列');
    });
    $('btnClear').addEventListener('click', function () {
      $('seqInput').value = ''; readInput();
    });
    $('btnFile').addEventListener('click', function () { $('fileInput').click(); });
    $('fileInput').addEventListener('change', function (e) {
      var f = e.target.files[0];
      if (!f) return;
      var fr = new FileReader();
      fr.onload = function () { $('seqInput').value = fr.result; readInput(); toast('已载入 ' + f.name); };
      fr.readAsText(f);
    });
    $('btnCopySeq').addEventListener('click', function () { copyText($('seqInput').value, '序列'); });
    $('btnDownloadSeq').addEventListener('click', function () {
      if (!state.seq) { toast('没有序列'); return; }
      var name = (state.name || 'sequence').replace(/\s+/g, '_');
      var content = '>' + name + ' | ' + state.seq.length + ' bp\n' + C.formatSeq(state.seq, 60) + '\n';
      download(name + '.fasta', content);
    });

    $('btnReport').addEventListener('click', function () {
      var md = buildReport();
      if (!md) { toast('请先输入序列'); return; }
      var base = (state.name || 'sequence').replace(/[^\w\u4e00-\u9fa5-]+/g, '_').slice(0, 40) || 'sequence';
      download(base + '-报告.md', md, 'text/markdown;charset=utf-8');
    });

    // 概览
    $('btnGcRedraw').addEventListener('click', drawGC);
    ['gcWin', 'gcStep'].forEach(function (id) { $(id).addEventListener('change', drawGC); });
    ['odVal', 'odDil', 'odKind'].forEach(function (id) { $(id).addEventListener('input', calcOD); });
    // 双链 MW 对碱基组成不敏感（A·T 对 617.41 vs G·C 对 618.39 Da），
    // 故这里用 AT 各半的等长序列估算即可，无需真实序列
    function approxDsMw(len) {
      var a = Math.floor(len / 2), t = Math.ceil(len / 2);
      return C.mwNucleic('A'.repeat(a) + 'T'.repeat(t), 'dsDNA');
    }
    $('btnCcFromNg').addEventListener('click', function () {
      var len = Math.max(1, parseInt($('ccLen').value, 10) || 1000);
      var ng = parseFloat($('ccNg').value);
      if (isNaN(ng)) { $('ccOut').textContent = '请填写 ng/µL'; return; }
      var mw = approxDsMw(len);
      var nm = C.concNgUlToNm(ng, mw);
      $('ccNm').value = Number(nm.toFixed(3));
      $('ccOut').innerHTML = len + ' bp 双链，MW ≈ ' + num(mw, 0) + ' Da：<b class="result-ok">' +
        num(ng, 2) + ' ng/µL = ' + nm.toFixed(2) + ' nM</b>（= ' + (nm / 1000).toFixed(4) + ' µM）';
    });
    $('btnCcFromNm').addEventListener('click', function () {
      var len = Math.max(1, parseInt($('ccLen').value, 10) || 1000);
      var nm = parseFloat($('ccNm').value);
      if (isNaN(nm)) { $('ccOut').textContent = '请填写 nM'; return; }
      var mw = approxDsMw(len);
      var ng = C.concNmToNgUl(nm, mw);
      $('ccNg').value = Number(ng.toFixed(3));
      $('ccOut').innerHTML = len + ' bp 双链，MW ≈ ' + num(mw, 0) + ' Da：<b class="result-ok">' +
        num(nm, 2) + ' nM = ' + ng.toFixed(3) + ' ng/µL</b>';
    });

    // 序列操作
    $('btnRevComp').addEventListener('click', function () { showOp(C.formatSeq(C.revComp(state.seq), 60)); });
    $('btnComp').addEventListener('click', function () {
      var m = { A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', N: 'N' };
      showOp(C.formatSeq(state.seq.replace(/[ATGCU]/g, function (c) { return m[c] || 'N'; }), 60));
    });
    $('btnRev').addEventListener('click', function () { showOp(C.formatSeq(state.seq.split('').reverse().join(''), 60)); });
    $('btnTranscribe').addEventListener('click', function () { showOp(C.formatSeq(C.transcribe(state.seq), 60)); });
    $('btnBackTranscribe').addEventListener('click', function () { showOp(C.formatSeq(C.reverseTranscribe(state.seq), 60)); });
    $('btnCopyOut').addEventListener('click', function () { copyText($('opOut').textContent.replace(/\n/g, ''), '结果'); });

    $('btnTranslate').addEventListener('click', function () {
      var v = $('frameSel').value;
      var s = state.seq;
      if (v.charAt(0) === 'r') s = C.revComp(s);
      var frame = parseInt(v.replace('r', ''), 10) || 0;
      var prot = C.translate(s, frame);
      $('protOut').innerHTML = colorAa(prot);
    });
    $('btnCopyProt').addEventListener('click', function () { copyText($('protOut').textContent.replace(/\n/g, ''), '蛋白序列'); });

    $('btnFindOrf').addEventListener('click', function () {
      // 同样不能用 || 兜底：用户填 0 会被换成 60
      var minAaRaw = parseInt($('orfMinAa').value, 10);
      var minAa = isNaN(minAaRaw) ? 60 : Math.max(1, minAaRaw);
      var orfs = getORFs(minAa);
      var html = '';
      if (!orfs.length) {
        html = '<tr><td colspan="7" style="color:var(--text-3)">未找到满足条件的 ORF，可降低最小长度或改用「不要求 ATG 起始」的思路检查序列方向</td></tr>';
      }
      orfs.slice(0, 60).forEach(function (o, i) {
        // 跨接缝的 ORF 起点大于终点（起点在序列末尾、绕回开头）
        var wrapTag = o.wrap
          ? ' <span class="tag warn" title="该 ORF 跨越序列首尾接缝，仅在环状模式下能检出">跨接缝</span>'
          : '';
        var posText = o.wrap ? (o.ntStart + ' → ' + o.ntEnd) : (o.ntStart + ' – ' + o.ntEnd);
        html += '<tr><td class="num">' + (i + 1) + '</td>' +
          '<td><span class="tag ' + (o.strand === '+' ? 'brand' : '') + '">' + o.strand + '</span></td>' +
          '<td class="num">' + o.frame + '</td>' +
          '<td class="num">' + posText + wrapTag + '</td>' +
          '<td class="num"><b>' + o.aaLength + '</b></td>' +
          '<td style="font-family:var(--mono);font-size:11.5px">' + esc(o.protein.slice(0, 40)) + (o.protein.length > 40 ? '…' : '') + '</td>' +
          '<td><button class="small link" data-orf="' + i + '">查看</button></td></tr>';
      });
      var tb = $('orfTbl').querySelector('tbody');
      tb.innerHTML = html;
      Array.prototype.forEach.call(tb.querySelectorAll('button[data-orf]'), function (b) {
        b.addEventListener('click', function () {
          var o = orfs[parseInt(b.getAttribute('data-orf'), 10)];
          $('orfProt').innerHTML = colorAa(o.protein);
          $('orfNt').textContent = C.formatSeq(o.nt, 60);
        });
      });
      if (orfs.length) {
        $('orfProt').innerHTML = colorAa(orfs[0].protein);
        $('orfNt').textContent = C.formatSeq(orfs[0].nt, 60);
      }
      toast(orfs.length ? '找到 ' + orfs.length + ' 个 ORF' : '未找到 ORF');
    });

    // 酶切
    $('enzSearch').addEventListener('input', renderEnzymeList);
    Array.prototype.forEach.call(document.querySelectorAll('[data-quick]'), function (b) {
      b.addEventListener('click', function () {
        var mode = b.getAttribute('data-quick');
        var names = [];
        if (mode === 'common') names = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SalI', 'PstI', 'KpnI', 'NotI'];
        else if (mode === 'mcs') names = ['EcoRI', 'SacI', 'KpnI', 'SmaI', 'BamHI', 'XbaI', 'SalI', 'PstI', 'SphI', 'HindIII'];
        else if (mode === 'single' || mode === 'zero') {
          var s = state.seq;
          if (!s) { toast('请先输入序列'); return; }
          E.list.forEach(function (en) {
            var c = C.digest(s, [en], state.circular).perEnzyme[0].count;
            if (mode === 'single' ? c === 1 : c === 0) names.push(en.name);
          });
          toast((mode === 'single' ? '单切点酶 ' : '零切点酶 ') + names.length + ' 种');
        }
        state.selectedEnzymes = names;
        $('enzSearch').value = '';
        renderEnzymeList();   // 选中集合被程序改动，需同步 checkbox 状态
        renderDigest();
      });
    });
    ['gelPct', 'gelPerEnzyme', 'gelLabel'].forEach(function (id) {
      $(id).addEventListener('change', function () { renderDigest(); });
    });
    $('btnGelPng').addEventListener('click', function () {
      var cv = $('gelCanvas');
      try {
        var a = document.createElement('a');
        a.href = cv.toDataURL('image/png');
        a.download = 'MolBench-gel.png';
        a.click();
        toast('已导出凝胶图');
      } catch (e) { toast('导出失败'); }
    });

    // 自动引物设计
    $('btnDesign').addEventListener('click', function () {
      var s = state.seq;
      if (!s || state.type === 'protein' || state.type === 'unknown') { toast('需要有效的核酸序列'); return; }
      var o = primerOpts();
      o.start = Math.max(0, (parseInt($('dzStart').value, 10) || 1) - 1);
      var endVal = parseInt($('dzEnd').value, 10);
      o.end = endVal ? Math.min(s.length, endVal) : s.length;
      o.minLen = parseInt($('dzMinLen').value, 10) || 18;
      o.maxLen = parseInt($('dzMaxLen').value, 10) || 27;
      o.minTm = parseFloat($('dzMinTm').value);
      o.maxTm = parseFloat($('dzMaxTm').value);
      o.minGc = parseFloat($('dzMinGc').value);
      o.maxGc = parseFloat($('dzMaxGc').value);
      o.maxProduct = parseInt($('dzMaxProd').value, 10) || 4000;
      o.topN = 10;
      if (isNaN(o.minTm)) o.minTm = 55;
      if (isNaN(o.maxTm)) o.maxTm = 65;
      if (isNaN(o.minGc)) o.minGc = 40;
      if (isNaN(o.maxGc)) o.maxGc = 60;

      var t0 = Date.now();
      var pairs = C.designPrimers(s, o);
      var ms = Date.now() - t0;

      var tb = $('dzTbl').querySelector('tbody');
      if (!pairs.length) {
        // 根据模板整体 GC 给出针对性建议，而不是一句笼统的"放宽条件"
        var tplGc = C.gcContent(s);
        var hint = '没有满足条件的引物对。';
        if (tplGc > 65) hint += '该模板 GC 含量高达 ' + num(tplGc, 1) + '%，建议把 GC% 上限放宽到 70–75%，' +
          '并在反应中加入 3–5% DMSO 或改用高 GC 专用缓冲液。';
        else if (tplGc < 35) hint += '该模板 GC 含量仅 ' + num(tplGc, 1) + '%，建议把 GC% 下限放宽到 30–35%，' +
          '并适当降低退火温度。';
        else hint += '可放宽 Tm（如 52–68 °C）或 GC%（35–65%）范围，或扩大引物长度区间。';
        tb.innerHTML = '<tr><td colspan="10" style="color:var(--text-3)">' + esc(hint) + '</td></tr>';
        $('dzOut').textContent = '耗时 ' + ms + ' ms。模板整体 GC 为 ' + num(tplGc, 1) + '%。';
        return;
      }
      var html = '';
      pairs.forEach(function (p, i) {
        var crossTag = p.touch3
          ? '<span class="tag bad" title="3′ 端参与交叉配对，最易形成引物二聚体">3′ ' + num(p.crossDg, 1) + '</span>'
          : (p.crossDg <= -5 ? '<span class="tag bad">' + num(p.crossDg, 1) + '</span>'
            : (p.crossDg <= -2.5 ? '<span class="tag warn">' + num(p.crossDg, 1) + '</span>'
              : '<span class="tag ok">' + num(p.crossDg, 1) + '</span>'));
        function gcTag(v) {
          // 模板局部 GC 过高／过低都会降低扩增效率，标出来提醒
          var cls = (v > 70 || v < 30) ? 'warn' : '';
          return cls ? '<span class="tag warn">' + num(v, 1) + '</span>' : num(v, 1);
        }
        html += '<tr><td class="num">' + (i + 1) + '</td>' +
          '<td style="font-family:var(--mono);font-size:11.5px">' + esc(p.fwd.seq) + '</td>' +
          '<td style="font-family:var(--mono);font-size:11.5px">' + esc(p.rev.seq) + '</td>' +
          '<td class="num"><b>' + p.size + '</b></td>' +
          '<td class="num">' + num(p.fwd.tm, 1) + ' / ' + num(p.rev.tm, 1) + '</td>' +
          '<td class="num">' + num(p.fwd.gc, 1) + ' / ' + num(p.rev.gc, 1) + '</td>' +
          '<td class="num">' + gcTag(p.fwd.localGc) + ' / ' + gcTag(p.rev.localGc) + '</td>' +
          '<td class="num">' + crossTag + '</td>' +
          '<td class="num"><b>' + num(p.score, 1) + '</b></td>' +
          '<td><button class="small link" data-dz="' + i + '">使用</button></td></tr>';
      });
      tb.innerHTML = html;
      Array.prototype.forEach.call(tb.querySelectorAll('button[data-dz]'), function (b) {
        b.addEventListener('click', function () {
          var p = pairs[parseInt(b.getAttribute('data-dz'), 10)];
          $('prFwd').value = p.fwd.seq;
          $('prRev').value = p.rev.seq;
          primerAutoFilled = true;
          renderPrimer();
          toast('已填入第 ' + (parseInt(b.getAttribute('data-dz'), 10) + 1) + ' 对引物，产物 ' + p.size + ' bp');
        });
      });
      $('dzOut').innerHTML = '共 <b>' + pairs.length + '</b> 对候选，耗时 ' + ms + ' ms。' +
        '评分综合考虑 Tm 接近 60 °C、GC 接近 50%、发夹与自身二聚体的 ΔG₃₇、模板内结合位点唯一性，' +
        '以及结合位点周边的模板局部 GC。「交叉 ΔG」为两引物交叉二聚体的自由能，' +
        '标注「3′」表示配对发生在 3′ 端（最危险，会被聚合酶延伸成二聚体产物）。';
    });

    // 引物
    $('btnPrimer').addEventListener('click', renderPrimer);
    ['prFwd', 'prRev', 'prConc', 'prNa', 'prMg', 'prDntp'].forEach(function (id) {
      $(id).addEventListener('change', renderPrimer);
    });
    $('btnPickFromTpl').addEventListener('click', function () {
      var s = state.seq;
      if (!s || s.length < 80) { toast('需要至少 80 bp 的模板序列'); return; }
      var f = s.substr(0, 22);
      var r = C.revComp(s.substr(-22));
      $('prFwd').value = f; $('prRev').value = r;
      renderPrimer();
      toast('已从模板两端各取 22 nt');
    });
    $('btnCopyPcr').addEventListener('click', function () { copyText($('pcrSeq').textContent.replace(/\n/g, ''), '产物序列'); });
    $('btnPcrToInput').addEventListener('click', function () {
      var t = $('pcrSeq').textContent.replace(/\n/g, '');
      if (!t || t === '—') { toast('没有产物'); return; }
      $('seqInput').value = '>PCR_product ' + t.length + 'bp\n' + C.formatSeq(t, 60);
      $('chkCircular').checked = false;
      readInput();
      toast('产物已载入，全长 ' + t.length + ' bp');
    });

    // 序列比对
    $('btnAlign').addEventListener('click', doAlign);
    $('btnAlnFillA').addEventListener('click', function () {
      if (!state.seq) { toast('当前没有序列'); return; }
      $('alnA').value = C.formatSeq(state.seq, 60);
      toast('已填入 A');
    });
    $('btnAlnFillB').addEventListener('click', function () {
      if (!state.seq) { toast('当前没有序列'); return; }
      $('alnB').value = C.formatSeq(state.seq, 60);
      toast('已填入 B');
    });
    $('btnAlnSwap').addEventListener('click', function () {
      var t = $('alnA').value; $('alnA').value = $('alnB').value; $('alnB').value = t;
      if ($('alnOut').textContent !== '—') doAlign();
    });
    $('btnAlnDemo').addEventListener('click', function () {
      $('alnA').value = C.formatSeq(DEMO_ALN_A, 60);
      $('alnB').value = C.formatSeq(DEMO_ALN_B, 60);
      doAlign();
      toast('已载入演示对（含缺失与错配）');
    });

    // 虚拟克隆
    fillEnzymeSelect($('clBackEnz'), true);
    fillEnzymeSelect($('clInsEnz'), true);
    $('clBackEnz').value = 'EcoRI';
    $('clInsEnz').value = 'EcoRI';
    refreshCloneFrags('back');
    refreshCloneFrags('ins');
    ['clBackSeq', 'clInsSeq'].forEach(function (id) {
      $(id).addEventListener('input', function () {
        refreshCloneFrags(id === 'clBackSeq' ? 'back' : 'ins');
      });
    });
    ['clBackEnz', 'clInsEnz'].forEach(function (id) {
      $(id).addEventListener('change', function () {
        refreshCloneFrags(id === 'clBackEnz' ? 'back' : 'ins');
      });
    });
    $('btnClBackFromSeq').addEventListener('click', function () {
      if (!state.seq) { toast('当前没有序列'); return; }
      $('clBackSeq').value = C.formatSeq(state.seq, 60);
      refreshCloneFrags('back');
      toast('已把当前序列填为载体');
    });
    $('btnClInsFromSeq').addEventListener('click', function () {
      if (!state.seq) { toast('当前没有序列'); return; }
      $('clInsSeq').value = C.formatSeq(state.seq, 60);
      refreshCloneFrags('ins');
      toast('已把当前序列填为插入片段');
    });
    $('btnClSwap').addEventListener('click', function () {
      var t = $('clBackSeq').value;
      $('clBackSeq').value = $('clInsSeq').value;
      $('clInsSeq').value = t;
      refreshCloneFrags('back');
      refreshCloneFrags('ins');
    });
    $('btnClone').addEventListener('click', doClone);
    $('btnClCopy').addEventListener('click', function () {
      if (!window.__cloneSeq) { toast('还没有重组序列'); return; }
      copyText(window.__cloneSeq, '重组序列');
    });
    $('btnClToInput').addEventListener('click', function () {
      if (!window.__cloneSeq) { toast('还没有重组序列'); return; }
      $('seqInput').value = '>recombinant\n' + C.formatSeq(window.__cloneSeq, 60);
      $('chkCircular').checked = true;   // 重组质粒按环状处理
      readInput();
      toast('已载入重组质粒（按环状处理）');
    });

    $('btnAlnMultiDemo').addEventListener('click', function () {
      var base = DEMO_ALN_A;
      var v1 = base.slice(0, 10) + base.slice(11);           // 缺失 1 nt
      var v2 = base.slice(0, 15) + 'A' + base.slice(15);     // 插入 1 nt
      var v3 = base.slice(0, 20) + 'T' + base.slice(21);     // 替换 1 nt
      $('alnA').value = [
        '>参考序列', base,
        '>克隆-1 缺失', v1,
        '>克隆-2 插入', v2,
        '>克隆-3 替换', v3
      ].join('\n');
      $('alnB').value = '';
      doAlign();
      toast('已载入 4 条演示序列');
    });

    // 键盘快捷键：Ctrl/Cmd + 数字切换面板，Ctrl/Cmd + Enter 执行当前面板主操作
    document.addEventListener('keydown', function (ev) {
      var mod = ev.ctrlKey || ev.metaKey;
      if (!mod) return;
      if (ev.key >= '1' && ev.key <= '9') {
        var btns = document.querySelectorAll('#tabs button');
        var idx = parseInt(ev.key, 10) - 1;
        if (btns[idx]) { btns[idx].click(); ev.preventDefault(); }
        return;
      }
      if (ev.key === 'Enter') {
        // 显式声明每个面板的主操作按钮，而不是靠 `.primary` 样式类去找
        // （早期实现依赖 class，结果「序列操作」面板根本没有 .primary 按钮，
        //   Ctrl+Enter 在那里面板静默失效）
        var MAIN_ACTION = {
          'panel-seqops': 'btnFindOrf',
          'panel-protein': 'btnProtAnalyze',
          'panel-primer': 'btnDesign',
          'panel-align': 'btnAlign',
          'panel-clone': 'btnClone'
        };
        var active = document.querySelector('.panel.active');
        var id = active && MAIN_ACTION[active.id];
        if (id && $(id)) { $(id).click(); ev.preventDefault(); }
      }
    });

    bindMapHover();
    bindCalc();
    bindProtein();
    bindDigestExtra();
    bindMut();
    bindLig();
    window.addEventListener('resize', function () {
      clearTimeout(window.__rz);
      window.__rz = setTimeout(function () { render(currentTab); }, 200);
    });
  }

  /**
   * 氨基酸着色 + 每 60 个换行。
   * 必须先分组再着色：若先插入 <span> 再按字符数换行，换行的 60 字符边界会
   * 落在标签中间（实测每 60 aa 就有 20+ 处把 `</span>` 切成 `</sp` + `an>`），
   * 导致 span 永不闭合、结构错乱。
   */
  function colorAa(prot) {
    var out = [];
    for (var i = 0; i < prot.length; i += 60) {
      out.push(prot.slice(i, i + 60).replace(/[A-Z*]/g, function (c) {
        return '<span class="aa-' + c + '">' + c + '</span>';
      }));
    }
    return out.join('\n');
  }

  /* =========================================================
   * 序列历史（localStorage，仅存本机）
   * ========================================================= */
  var HIST_KEY = 'molbench.history.v1';
  var HIST_MAX = 10;
  // localStorage 常见配额约 5 MB，且 file:// 下可能更小。
  // 单条放宽到 100 kb 时，10 条就是 1 MB 量级，容易写满导致静默失败；
  // 这里收紧到 20 kb / 条，并对整体字节数设上限，超出时丢弃最旧的记录。
  var HIST_MAX_CHARS = 20000;
  var HIST_MAX_BYTES = 800000;

  function histAvailable() {
    try {
      localStorage.setItem('__molbench_probe__', '1');
      localStorage.removeItem('__molbench_probe__');
      return true;
    } catch (e) { return false; }
  }

  function histLoad() {
    try {
      var raw = localStorage.getItem(HIST_KEY);
      var list = raw ? JSON.parse(raw) : [];
      return Object.prototype.toString.call(list) === '[object Array]' ? list : [];
    } catch (e) { return []; }
  }

  function histAdd(name, seq, circular) {
    if (!seq || seq.length > HIST_MAX_CHARS || !histAvailable()) return;
    var list = histLoad().filter(function (it) { return it && it.seq !== seq; });
    list.unshift({ name: name, seq: seq, circular: !!circular, ts: Date.now() });
    list = list.slice(0, HIST_MAX);
    // 先按总字节数裁剪，再写入；配额仍不足时继续逐条丢弃最旧的
    while (list.length) {
      var json = JSON.stringify(list);
      if (json.length > HIST_MAX_BYTES) { list.pop(); continue; }
      try { localStorage.setItem(HIST_KEY, json); break; }
      catch (e) { list.pop(); }
    }
  }

  function histClear() {
    try { localStorage.removeItem(HIST_KEY); } catch (e) { /* 忽略：存储不可用时静默降级 */ }
  }

  function timeAgo(ts) {
    var d = Date.now() - ts;
    if (d < 60000) return '刚刚';
    if (d < 3600000) return Math.floor(d / 60000) + ' 分钟前';
    if (d < 86400000) return Math.floor(d / 3600000) + ' 小时前';
    return Math.floor(d / 86400000) + ' 天前';
  }

  function renderHistory() {
    var sel = $('histSel');
    if (!sel) return;
    var list = histLoad();
    sel.innerHTML = '';
    var head = document.createElement('option');
    head.value = '';
    head.textContent = list.length ? '— 历史记录（' + list.length + '）—' : '— 历史记录（空）—';
    sel.appendChild(head);
    list.forEach(function (it, i) {
      var o = document.createElement('option');
      o.value = String(i);
      var label = (it.name || '(未命名)').slice(0, 28);
      o.textContent = label + ' · ' + it.seq.length + ' nt · ' + timeAgo(it.ts);
      sel.appendChild(o);
    });
    sel.disabled = list.length === 0;
  }

  var histTimer = null;
  function histSchedule() {
    clearTimeout(histTimer);
    // 停止输入 2 秒后才入库，避免打字过程中反复写入
    histTimer = setTimeout(function () {
      if (state.seq && state.seq.length >= 30 && state.type !== 'unknown') {
        histAdd(state.name, state.seq, state.circular);
        renderHistory();
      }
    }, 2000);
  }

  /* =========================================================
   * 分析报告导出（Markdown）
   * ========================================================= */
  function typeLabel() {
    return { dna: 'DNA', rna: 'RNA', protein: '蛋白质', unknown: '未知' }[state.type] || '未知';
  }

  function buildReport() {
    var s = state.seq;
    if (!s) return null;
    var L = [];
    var isNt = (state.type === 'dna' || state.type === 'rna');

    // 章节编号由计数器生成：FEATURES 与比对章节是条件出现的，
    // 写死「一、二、三」在增删章节时必然错位。
    var secNo = 0;
    var SEC_CN = ['一', '二', '三', '四', '五', '六', '七', '八', '九', '十', '十一', '十二', '十三', '十四'];
    function nextSec(title) {
      secNo++;
      L.push('## ' + (SEC_CN[secNo - 1] || secNo) + '、' + title);
      L.push('');
    }
    var stamp = new Date();
    function pad(v) { return (v < 10 ? '0' : '') + v; }
    var timeStr = stamp.getFullYear() + '-' + pad(stamp.getMonth() + 1) + '-' + pad(stamp.getDate()) +
      ' ' + pad(stamp.getHours()) + ':' + pad(stamp.getMinutes());

    L.push('# 序列分析报告 · MolBench');
    L.push('');
    L.push('| 项目 | 值 |');
    L.push('| --- | --- |');
    L.push('| 生成时间 | ' + timeStr + ' |');
    L.push('| 序列名称 | ' + mdCell(state.name || '(未命名)') + ' |');
    L.push('| 序列类型 | ' + typeLabel() + ' |');
    L.push('| 序列长度 | ' + s.length + ' ' + (state.type === 'protein' ? 'aa' : 'nt') + ' |');
    if (isNt) {
      L.push('| GC 含量 | ' + num(C.gcContent(s), 2) + '% |');
      L.push('| 拓扑结构 | ' + (state.circular ? '环状（质粒 / 载体）' : '线性') + ' |');
    }
    L.push('');

    // 一、分子量
    nextSec('分子量与浓度');
    if (isNt) {
      L.push('| 项目 | 数值 |');
      L.push('| --- | --- |');
      L.push('| 单链 DNA（5′-OH） | ' + fmtDa(C.mwNucleic(s, 'ssDNA')) + ' |');
      L.push('| 双链 DNA | ' + fmtDa(C.mwNucleic(s, 'dsDNA')) + ' |');
      L.push('| 单链 RNA | ' + fmtDa(C.mwNucleic(s, 'ssRNA')) + ' |');
      var mw = C.mwNucleic(s, 'dsDNA');
      L.push('| 50 ng/µL 对应 | ' + num(C.concNgUlToNm(50, mw), 2) + ' nM |');
    } else if (state.type === 'protein') {
      L.push('蛋白质分子量：**' + fmtDa(C.mwProtein(s)) + '**');
    }
    L.push('');

    // 二、ORF
    if (isNt) {
      nextSec('开放阅读框（≥60 aa）');
      var orfs = getORFs(60);
      if (!orfs.length) {
        L.push('未检出满足条件的 ORF。');
      } else {
        L.push('| # | 链 | 框 | 区间 (nt) | 长度 (aa) | 起始序列 |');
        L.push('| --- | --- | --- | --- | --- | --- |');
        orfs.slice(0, 15).forEach(function (o, i) {
          L.push('| ' + (i + 1) + ' | ' + o.strand + ' | ' + o.frame + ' | ' +
            o.ntStart + '–' + o.ntEnd + ' | ' + o.aaLength + ' | `' + o.protein.slice(0, 20) + '…` |');
        });
        if (orfs.length > 15) L.push('');
        if (orfs.length > 15) L.push('（共 ' + orfs.length + ' 个，仅列出前 15 个）');
      }
      L.push('');
    }

    // 序列注释（GenBank FEATURES）
    var feats = state.features || [];
    if (feats.length) {
      nextSec('序列注释（GenBank FEATURES）');
      L.push('| 类型 | 位置 | 区间 | 链 | 基因 / 产物 |');
      L.push('| --- | --- | --- | --- | --- |');
      feats.forEach(function (f) {
        var label = f.qualifiers.gene || f.qualifiers.product || f.qualifiers.locus_tag || '—';
        if (f.qualifiers.gene && f.qualifiers.product) label = f.qualifiers.gene + ' — ' + f.qualifiers.product;
        L.push('| ' + mdCell(f.type) + ' | `' + mdCell(f.location) + '` | ' +
          (f.start ? f.start + '–' + f.end : '—') + ' | ' + f.strand + ' | ' + mdCell(label) + ' |');
      });
      L.push('');
    }

    // 酶切
    if (isNt) {
      nextSec('限制性酶切分析');
      var sel = state.selectedEnzymes.map(function (n) { return E.byName(n); }).filter(Boolean);
      if (!sel.length) {
        L.push('未选择限制性内切酶。');
      } else {
        var d = C.digest(s, sel, state.circular);
        L.push('**使用酶**：' + sel.map(function (e) { return e.name + '（`' + e.site + '`）'; }).join('、'));
        L.push('');
        L.push('**反应拓扑**：' + (state.circular ? '环状' : '线性') + '，共 ' + d.fragments.length + ' 个片段');
        L.push('');
        L.push('| # | 酶 | 识别序列 | 切点位置 | 备注 |');
        L.push('| --- | --- | --- | --- | --- |');
        var rows = 0;
        d.perEnzyme.forEach(function (pe) {
          pe.hits.forEach(function (h) {
            rows++;
            L.push('| ' + rows + ' | ' + pe.enzyme.name + ' | `' + pe.enzyme.site + '` | ' +
              h.cut + '/' + (h.cut + 1) + ' | ' + (h.wrap ? '跨接缝' : '') + ' |');
          });
        });
        if (!rows) L.push('| — | — | — | — | 所选酶无切点 |');
        L.push('');
        L.push('**片段长度**：' + d.fragments.join(' bp、') + ' bp');
        L.push('');
        L.push('> 合计 ' + d.fragments.reduce(function (a, b) { return a + b; }, 0) + ' bp。');
      }
      L.push('');
    }

    // 四、引物与 PCR
    nextSec('引物与 PCR');
    var fwd = ($('prFwd').value || '').toUpperCase().replace(/\s/g, '');
    var rev = ($('prRev').value || '').toUpperCase().replace(/\s/g, '');
    if (fwd || rev) {
      var opts = primerOpts();
      L.push('| 引物 | 序列 (5′→3′) | 长度 | GC% | Tm 最近邻 | 3′ 端 |');
      L.push('| --- | --- | --- | --- | --- | --- |');
      [{ n: '正向', v: fwd }, { n: '反向', v: rev }].forEach(function (p) {
        if (!p.v) return;
        L.push('| ' + p.n + ' | `' + p.v + '` | ' + p.v.length + ' | ' + num(C.gcContent(p.v), 1) +
          ' | ' + num(C.tmNearestNeighbor(p.v, opts), 1) + ' °C | `' + C.gcClamp(p.v).tail + '` |');
      });
      L.push('');
      if (fwd && rev && isNt) {
        var prs = C.pcrProducts(s, fwd, rev, { circular: state.circular });
        if (prs.length === 1) {
          L.push('**产物大小**：' + prs[0].size + ' bp（' + prs[0].start + '–' + prs[0].end + '）');
        } else if (prs.length > 1) {
          L.push('**模板存在多个引物结合位点**，检出 ' + prs.length + ' 种可能产物（按长度升序，越短越易被优先扩增）：');
          L.push('');
          L.push('| # | 产物 (bp) | 区间 | F 位点 | R 位点 |');
          L.push('| --- | --- | --- | --- | --- |');
          prs.forEach(function (p, i) {
            L.push('| ' + (i + 1) + ' | ' + p.size + ' | ' + p.start + '–' + p.end + ' | #' + p.fwdIndex + ' | #' + p.revIndex + ' |');
          });
        } else {
          L.push('未能在模板上定位有效的引物结合位点组合。');
        }
        if (prs.length) {
          var ta = Math.min(C.tmNearestNeighbor(fwd, opts), C.tmNearestNeighbor(rev, opts)) - 3;
          L.push('');
          L.push('**建议退火温度**：' + ta.toFixed(1) + ' °C；延伸时间按 1 kb/min 估算为 ' +
            Math.max(15, Math.ceil(prs[0].size / 1000) * 60) + ' s。');
        }
      }
    } else {
      L.push('未输入引物。');
    }
    L.push('');

    // 双序列比对（用户填过才输出，避免报告里出现空白章节）
    var alnA = ($('alnA').value || '').trim();
    var alnB = ($('alnB').value || '').trim();
    if (alnA && alnB) {
      function intOf(id, dflt) {
        var v = parseInt($(id).value, 10);
        return isNaN(v) ? dflt : v;
      }
      var ar = C.alignPair(alnA, alnB, {
        match: intOf('alnMatch', 2),
        mismatch: intOf('alnMis', -1),
        gap: intOf('alnGap', -2)
      });
      nextSec('双序列比对（Needleman-Wunsch）');
      if (ar && ar.tooBig) {
        L.push('序列过长（需 ' + ar.cells + ' 个矩阵单元），未在报告中展开比对结果。');
      } else if (ar) {
        L.push('| 项目 | 值 |');
        L.push('| --- | --- |');
        L.push('| 一致度 | ' + num(ar.identity, 2) + '% |');
        L.push('| 比对得分 | ' + ar.score + ' |');
        L.push('| 一致 / 错配 / gap | ' + ar.match + ' / ' + (ar.length - ar.match - ar.gaps) + ' / ' + ar.gaps + ' |');
        L.push('| 比对长度 (A : B) | ' + ar.length + '（' + ar.lenA + ' : ' + ar.lenB + '） |');
      } else {
        L.push('比对输入为空或只含非法字符。');
      }
      L.push('');
    }

    // 蛋白理化性质与密码子优化：蛋白面板填过内容才输出，避免报告出现空章节
    var repProt = readProtSeq().replace(/\*/g, '');
    var pp = repProt ? C.proteinProperties(repProt) : null;
    if (pp) {
      nextSec('蛋白理化性质');
      L.push('| 指标 | 值 |');
      L.push('| --- | --- |');
      L.push('| 长度 | ' + pp.length + ' aa |');
      L.push('| 分子量 | ' + num(pp.mw, 1) + ' Da |');
      L.push('| 等电点 pI | ' + num(pp.pI, 2) + ' |');
      L.push('| pH 7 净电荷 | ' + (pp.charge7 >= 0 ? '+' : '') + num(pp.charge7, 2) + ' |');
      L.push('| GRAVY 疏水性 | ' + num(pp.gravy, 3) + ' |');
      L.push('| 脂肪族指数 | ' + num(pp.aliphatic, 2) + ' |');
      L.push('| 不稳定指数 | ' + num(pp.instability, 2) + (pp.unstable ? '（>40，预测不稳定）' : '（≤40，预测稳定）') + ' |');
      L.push('| 消光系数 ε₂₈₀ | ' + num(pp.extinctionOx, 0) + ' M⁻¹cm⁻¹（全氧化）/ ' + num(pp.extinctionRed, 0) + '（还原） |');
      L.push('| A₂₈₀ 0.1% | ' + num(pp.abs01Ox, 3) + ' |');
      L.push('');
      L.push('> pI 用 EMBOSS/Bjellqvist pKa + Henderson-Hasselbalch 净电荷二分求解；GRAVY 用 Kyte-Doolittle；');
      L.push('> 不稳定指数用 Guruprasad (1990) DIWV 二肽权重；消光系数用 Pace (1995)。');
      L.push('');

      var rt = C.reverseTranslate(repProt);
      nextSec('密码子优化（E. coli K12 偏好密码子）');
      L.push('- DNA 长度：' + rt.dna.length + ' bp；GC 含量：' + num(C.gcContent(rt.dna), 2) + '%');
      L.push('- 回译校验：' + (rt.backTranslateOk ? '与原蛋白序列完全一致' : '**不一致，请检查输入**') +
        (rt.skipped ? '；跳过非常规残基 ' + rt.skipped + ' 个' : ''));
      var optHits = [];
      E.list.forEach(function (en) {
        var hh = C.findSites(rt.dna, en, false);
        if (hh.length) optHits.push(en.name + '×' + hh.length);
      });
      L.push('- 产物内常用酶切位点：' + (optHits.length ? optHits.join('、') + '（如需亚克隆建议用同义密码子规避）' : '51 种常用酶均无切点'));
      L.push('');
      L.push('`' + rt.dna + '`');
      L.push('');
    }

    // 定点突变引物：填过突变位置才输出
    if (isNt && ($('mutPos').value || '').trim()) {
      var mres = C.designMutPrimers(s, parseInt($('mutPos').value, 10), $('mutNew').value, {
        oldLen: parseInt($('mutOld').value, 10) || 0,
        targetTm: parseFloat($('mutTm').value) || 62,
        circular: state.circular,
        dnaConc: parseFloat($('prConc').value) || 500,
        naConc: parseFloat($('prNa').value) || 50
      });
      nextSec('定点突变引物（Q5® SDM 型）');
      if (!mres.ok) {
        L.push('未能设计：' + mres.reason);
      } else {
        L.push('| 项目 | 内容 |');
        L.push('| --- | --- |');
        L.push('| 突变 | ' + mres.mutDesc + ' |');
        L.push('| 正向引物 (5′→3′) | `' + mres.fwd + '`（红字突变部 `' + (mres.fwdTail || '无') + '`，退火区 Tm ' + num(mres.tmF, 1) + ' °C） |');
        L.push('| 反向引物 (5′→3′) | `' + mres.rev + '`（退火区 Tm ' + num(mres.tmR, 1) + ' °C） |');
        L.push('| 突变后质粒大小 | ' + mres.productLen + ' bp |');
        L.push('');
        L.push('> 流程：PCR 扩增全长质粒 → KLD（激酶-连接酶-DpnI）处理 → 转化。');
        mres.warnings.forEach(function (w) { L.push('> ⚠ ' + w); });
      }
      L.push('');
    }

    // Motif 搜索：填过模式才输出
    var repPat = ($('motifPat').value || '').toUpperCase().replace(/\s/g, '');
    if (isNt && repPat && /^[ACGTURYSWKMBDHVN]+$/.test(repPat)) {
      repPat = repPat.replace(/U/g, 'T');
      var mhits = C.findMotifs(s, repPat, {
        circular: state.circular,
        bothStrands: $('motifBoth').checked
      });
      nextSec('Motif 搜索（' + repPat + '）');
      L.push('共 ' + mhits.length + ' 处命中' + (mhits.truncated ? '（已达上限被截断）' : '') +
        '；搜索范围：' + (state.circular ? '环状' : '线性') + ($('motifBoth').checked ? '、双链' : '、仅正链') + '。');
      L.push('');
      if (mhits.length) {
        L.push('| # | 位置 | 链 | 命中序列 | 备注 |');
        L.push('| --- | --- | --- | --- | --- |');
        mhits.slice(0, 50).forEach(function (hh, i) {
          L.push('| ' + (i + 1) + ' | ' + hh.pos + ' | ' + (hh.strand === '+' ? '正链' : '负链') + ' | `' + hh.match + '` | ' + (hh.wrap ? '跨接缝' : '') + ' |');
        });
        if (mhits.length > 50) L.push('');
        if (mhits.length > 50) L.push('> 仅列出前 50 条，完整结果请在界面中查看。');
      }
      L.push('');
    }

    L.push('---');
    L.push('');
    L.push('> 本报告由 MolBench（离线分子生物实验台）生成。Tm 采用 SantaLucia (1998) 最近邻参数 + Owczarzy (2004) 盐校正；');
    L.push('> 虚拟凝胶迁移距离为经验模型。结果供实验设计参考，关键实验请以实测为准。');
    return L.join('\n');
  }

  function download(filename, content, mime) {
    var blob = new Blob([content], { type: mime || 'text/plain;charset=utf-8' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a); a.click();
    setTimeout(function () { document.body.removeChild(a); URL.revokeObjectURL(a.href); }, 200);
    toast('已导出 ' + filename);
  }

  /* =========================================================
   * 启动
   * ========================================================= */
  function init() {
    bind();
    // 酶列表原先由 renderDigest() 顺带渲染；解耦后必须在这里显式初始化，
    // 否则首次进入酶切面板时列表是空的。
    renderEnzymeList();
    if (EX.length) {
      var e = EX[0];
      $('seqInput').value = '>' + e.id + ' ' + e.name + '\n' + C.formatSeq(e.seq, 60);
      $('chkCircular').checked = !!e.circular;
    }
    readInput();

    renderHistory();

    // 支持 URL hash 直达面板，例如 index.html#digest
    var hash = (location.hash || '').replace('#', '');
    if (hash) {
      var btn = document.querySelector('#tabs button[data-panel="' + hash + '"]');
      if (btn) btn.click();
    }
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
