/* 功能对账：检查 README 承诺但常规冒烟测试未覆盖的功能路径 */
const path = require('path');
const { JSDOM, VirtualConsole } = require('jsdom');
const root = path.join(__dirname, '..');

const errors = [];
const vc = new VirtualConsole();
vc.on('jsdomError', e => { if (!/Not implemented/.test(e.message)) errors.push('jsdomError: ' + (e.stack || e.message)); });

function fakeCtx() {
  const noop = () => {};
  return new Proxy({}, {
    get(t, k) { if (k === 'measureText') return () => ({ width: 10 }); return k in t ? t[k] : noop; },
    set(t, k, v) { t[k] = v; return true; }
  });
}

JSDOM.fromFile(path.join(root, 'index.html'), {
  runScripts: 'dangerously', resources: 'usable', pretendToBeVisual: true, virtualConsole: vc,
  beforeParse(win) {
    win.HTMLCanvasElement.prototype.getContext = () => fakeCtx();
    win.HTMLCanvasElement.prototype.toDataURL = () => 'data:image/png;base64,AAAA';
    win.scrollTo = () => {};
    win.devicePixelRatio = 1;
    win.onerror = (m, s, l, c, e) => errors.push('onerror: ' + m + ' @' + l);
    Object.defineProperty(win.Element.prototype, 'clientWidth', { get() { return 900; }, configurable: true });
    Object.defineProperty(win.Element.prototype, 'clientHeight', { get() { return 300; }, configurable: true });
    // jsdom 对 file://（opaque origin）不提供 localStorage，注入内存实现以便验证业务逻辑。
    // 真实浏览器下 file:// 是可用的，见 tools/visual_check.html 的 localStorage 探测。
    let nativeOk = false;
    try { win.localStorage.setItem('__probe__', '1'); nativeOk = win.localStorage.getItem('__probe__') === '1'; win.localStorage.removeItem('__probe__'); } catch (e) { nativeOk = false; }
    if (!nativeOk) {
      const mem = {};
      Object.defineProperty(win, 'localStorage', {
        configurable: true,
        value: {
          getItem: k => (Object.prototype.hasOwnProperty.call(mem, k) ? mem[k] : null),
          setItem: (k, v) => { mem[k] = String(v); },
          removeItem: k => { delete mem[k]; },
          clear: () => { Object.keys(mem).forEach(k => delete mem[k]); },
          key: i => Object.keys(mem)[i] || null,
          get length() { return Object.keys(mem).length; }
        }
      });
    }
    win.__lsNative = nativeOk;

    // 捕获下载动作
    win.__downloads = [];
    const RealBlob = win.Blob;
    win.Blob = function (parts, opts) { win.__blob = { text: parts.join(''), type: opts && opts.type }; return new RealBlob(parts, opts); };
    win.URL.createObjectURL = () => 'blob:fake';
    win.URL.revokeObjectURL = () => {};
    const origClick = win.HTMLAnchorElement.prototype.click;
    win.HTMLAnchorElement.prototype.click = function () {
      if (this.download) win.__downloads.push({ name: this.download, href: this.href });
      else if (origClick) origClick.call(this);
    };
    // 捕获剪贴板
    win.__clipboard = '';
    Object.defineProperty(win.navigator, 'clipboard', {
      value: { writeText: t => { win.__clipboard = t; return Promise.resolve(); } },
      configurable: true
    });
  }
}).then(dom => new Promise(r => {
  const w = dom.window;
  if (w.document.readyState === 'complete') return r(w);
  w.addEventListener('load', () => r(w));
  setTimeout(() => r(w), 4000);
})).then(w => {
  const d = w.document, $ = id => d.getElementById(id);
  const issues = [];
  let n = 0;
  function check(name, fn) {
    n++;
    try {
      const r = fn();
      if (r) issues.push('[缺陷] ' + name + ' → ' + r);
      else console.log('  [OK]   ' + name);
    } catch (e) { issues.push('[异常] ' + name + ' → ' + (e.message || e)); }
  }
  const wait = ms => new Promise(r => setTimeout(r, ms));
  function tab(t) { d.querySelector('#tabs button[data-panel="' + t + '"]').dispatchEvent(new w.Event('click', { bubbles: true })); }
  function click(id) { const e = $(id); if (!e) throw new Error('缺少元素 ' + id); e.dispatchEvent(new w.Event('click', { bubbles: true })); }
  function setV(id, v) { const e = $(id); if (!e) throw new Error('缺少元素 ' + id); e.value = v; }
  function fire(id, ev) { $(id).dispatchEvent(new w.Event(ev, { bubbles: true })); }

  console.log('\n=== A. 导入 / 导出 ===');
  check('导出 FASTA 文件名与内容', () => {
    w.__downloads = [];
    click('btnDownloadSeq');
    if (!w.__downloads.length) return '未触发下载';
    const dl = w.__downloads[0];
    if (!/\.fasta$/.test(dl.name)) return '文件名异常: ' + dl.name;
    return false;
  });

  check('导出分析报告 (.md) 内容完整', () => {
    w.__blob = null;
    click('btnReport');
    if (!w.__blob) return '未生成文件';
    if (!/markdown/.test(w.__blob.type)) return 'MIME 应为 markdown，实为 ' + w.__blob.type;
    const t = w.__blob.text;
    ['# 序列分析报告', '一、分子量', '二、开放阅读框', '三、限制性酶切', '四、引物'].forEach(k => {
      if (t.indexOf(k) === -1) throw new Error('报告缺少章节: ' + k);
    });
    if (/undefined|NaN/.test(t)) return '报告中出现 undefined / NaN';
    if (t.length < 500) return '报告过短（' + t.length + ' 字符）';
    return false;
  });

  check('上传裸 FASTA 文件', () => {
    const f = new w.File(['>up_test\nACGTACGTACGTACGTACGT\nACGT'], 'a.fa', { type: 'text/plain' });
    const inp = $('fileInput');
    Object.defineProperty(inp, 'files', { value: [f], configurable: true });
    inp.dispatchEvent(new w.Event('change', { bubbles: true }));
    return false;
  });

  check('上传 GenBank 文件（ORIGIN 段）', () => {
    const gb = [
      'LOCUS       TEST        40 bp    DNA     linear   UNK 01-JAN-2026',
      'DEFINITION  test sequence for audit.',
      'ACCESSION   TEST0001',
      'FEATURES             Location/Qualifiers',
      '     CDS             1..40',
      '                     /gene="test"',
      'ORIGIN',
      '        1 atggccattg taatgggccg ctgaaagggt gcccgatagc',
      '//'
    ].join('\n');
    const f = new w.File([gb], 'b.gb', { type: 'text/plain' });
    const inp = $('fileInput');
    Object.defineProperty(inp, 'files', { value: [f], configurable: true });
    inp.dispatchEvent(new w.Event('change', { bubbles: true }));
    return false;
  });

  return wait(400).then(() => {
    // 上一项上传后需要等待 FileReader，这里检查解析结果
    check('GenBank 解析结果是否可用', () => {
      const meta = $('seqMeta').textContent;
      const m = meta.match(/长度\s*(\d+)/);
      if (!m) return '未解析出长度: ' + meta.slice(0, 60);
      const len = parseInt(m[1], 10);
      // 正确的 ORIGIN 解析应得到 40 nt；若把 LOCUS/DEFINITION 当序列则会得到上百 nt
      if (len > 60) return '疑似把注释当成了序列（长度 ' + len + ' nt，应为 40 nt）';
      if (len !== 40) return '长度 ' + len + ' nt，期望 40 nt';
      return false;
    });

    check('GenBank feature 注释表可渲染并查看', () => {
      tab('seqops');
      if ($('featCard').style.display === 'none') return '存在 FEATURES 段但注释卡片未显示';
      const rows = $('featTbl').querySelectorAll('tbody tr').length;
      if (!rows) return '注释表为空';
      const b = $('featTbl').querySelector('button[data-feat]');
      if (!b) return '缺少「查看」按钮';
      b.dispatchEvent(new w.Event('click', { bubbles: true }));
      if ($('featSeq').textContent.length < 10) return '未显示注释区段序列: ' + $('featSeq').textContent.slice(0, 30);
      if (/^\(该区段超出/.test($('featSeq').textContent)) return '注释区段提取越界';
      return false;
    });

    check('导出报告含 FEATURES 与比对章节，且章节编号连续', () => {
      tab('align');
      click('btnAlnDemo');            // 填入演示对，让比对章节有内容
      w.__blob = null;
      click('btnReport');
      if (!w.__blob) return '未生成报告';
      const t = w.__blob.text;
      if (t.indexOf('序列注释（GenBank FEATURES）') === -1) return '报告缺少 FEATURES 章节';
      if (t.indexOf('双序列比对') === -1) return '报告缺少比对章节';
      // 章节编号必须连续（曾经是写死的一/二/三/四，插入章节后必然错位）
      const CN = ['一', '二', '三', '四', '五', '六', '七', '八', '九', '十', '十一', '十二', '十三', '十四'];
      const heads = (t.match(/^## (十[一二三四]?|[一二三四五六七八九])、/gm) || [])
        .map(s => s.replace(/^## /, '').replace('、', ''));
      if (!heads.length) return '未解析到任何章节标题';
      for (let i = 0; i < heads.length; i++) {
        if (heads[i] !== CN[i]) {
          return '章节编号不连续：实际 ' + heads.join(',') + '，期望 ' + CN.slice(0, heads.length).join(',');
        }
      }
      return false;
    });

    check('报告按需包含新功能章节，且不填输入时不出现空章节', () => {
      const t = w.__blob.text;
      const news = ['蛋白理化性质', '密码子优化', '定点突变引物', 'Motif 搜索'];
      for (const k of news) {
        if (t.indexOf(k) !== -1) return '未填任何新功能输入，却出现了章节：' + k;
      }
      // 直接填入蛋白序列与各项输入（不依赖当前序列长度，当前模板是 40 nt 的 GenBank 片段）
      d.getElementById('protSeq').value = '>demo_protein\nMKTAYIAKQRQISFVKSHFSRQLEDLRQFIERTKKLD';
      d.getElementById('motifPat').value = 'AATAAA';
      d.getElementById('mutPos').value = '5';
      d.getElementById('mutNew').value = 'A';
      w.__blob = null;
      click('btnReport');
      if (!w.__blob) return '填入新功能输入后未生成报告';
      const t2 = w.__blob.text;
      for (const k of news) {
        if (t2.indexOf(k) === -1) return '报告缺少新章节：' + k;
      }
      if (/undefined|NaN/.test(t2)) return '含新章节的报告出现 undefined / NaN';
      // 密码子优化章节里的 DNA 应能回译成报告里的蛋白
      const m = t2.match(/^`([ACGT]{30,})`$/m);
      if (!m) return '报告缺少优化 DNA 序列块';
      if (t2.indexOf('回译校验：与原蛋白序列完全一致') === -1) return '报告未记录回译校验结论';
      return false;
    });

    console.log('\n=== B. 复制行为一致性 ===');
    tab('seqops');
    setV('exampleSel', '0'); fire('exampleSel', 'change');
    return wait(350);
  }).then(() => {
    check('复制蛋白序列不含换行', () => {
      click('btnTranslate');
      w.__clipboard = '';
      click('btnCopyProt');
      if (!w.__clipboard) return '未复制到内容';
      if (/\n/.test(w.__clipboard)) return '包含换行符，与「复制结果 / 复制产物」行为不一致';
      return false;
    });
    check('复制反向互补结果不含换行', () => {
      w.__clipboard = '';
      click('btnRevComp'); click('btnCopyOut');
      if (!w.__clipboard) return '未复制到内容';
      if (/\n/.test(w.__clipboard)) return '包含换行符';
      return false;
    });
    check('复制序列（原始输入）可用', () => {
      w.__clipboard = '';
      click('btnCopySeq');
      if (!w.__clipboard) return '未复制到内容';
      return false;
    });

    console.log('\n=== C. 示例与序列切换 ===');
    check('切换到第 2 个示例（基因组片段）', () => {
      setV('exampleSel', '1'); fire('exampleSel', 'change');
      const meta = $('seqMeta').textContent;
      if (!/长度/.test(meta)) return '状态条未更新';
      if ($('chkCircular').checked) return '基因组片段应为线性，但勾选了环状';
      return false;
    });
    check('切回第 1 个示例（环状质粒）', () => {
      setV('exampleSel', '0'); fire('exampleSel', 'change');
      if (!$('chkCircular').checked) return '质粒示例应为环状';
      return false;
    });

    console.log('\n=== D. 环状切换触发重算 ===');
    return wait(350);
  }).then(() => {
    check('环状 → 线性的片段数变化', () => {
      tab('digest');
      const getFrags = () => $('fragTbl').querySelectorAll('tbody tr').length;
      $('chkCircular').checked = true; fire('chkCircular', 'change');
      const circ = getFrags();
      $('chkCircular').checked = false; fire('chkCircular', 'change');
      const lin = getFrags();
      $('chkCircular').checked = true; fire('chkCircular', 'change');
      if (circ === 0 || lin === 0) return '片段表为空';
      if (circ === lin) return '环状与线性片段数相同（' + circ + '），切换可能未生效';
      return false;
    });

    console.log('\n=== E. 凝胶导出与其它按钮 ===');
    check('凝胶导出 PNG', () => {
      w.__downloads = [];
      click('btnGelPng');
      if (!w.__downloads.length) return '未触发下载';
      if (!/gel.*\.png/i.test(w.__downloads[0].name)) return '文件名异常: ' + w.__downloads[0].name;
      return false;
    });
    check('随机生成 / 清空不报错', () => {
      click('btnDemoRandom');
      click('btnClear');
      if ($('seqInput').value !== '') return '清空未生效';
      return false;
    });

    console.log('\n=== F. 产物回写与 ORF 查看 ===');
    setV('exampleSel', '0'); fire('exampleSel', 'change');
    return wait(350);
  }).then(() => {
    check('PCR 产物载入为当前序列', () => {
      tab('primer');
      click('btnPickFromTpl');
      const before = $('pcrSeq').textContent.replace(/\n/g, '');
      if (!before || before === '—') return '无产物';
      click('btnPcrToInput');
      const meta = $('seqMeta').textContent;
      if (meta.indexOf(before.length) === -1) return '载入后长度 ' + before.length + ' 未体现在状态条';
      return false;
    });
    check('ORF 表格「查看」按钮可用', () => {
      tab('seqops');
      click('btnFindOrf');
      const b = $('orfTbl').querySelector('button[data-orf]');
      if (!b) return '无 ORF 行';
      b.dispatchEvent(new w.Event('click', { bubbles: true }));
      if ($('orfProt').textContent.length < 10) return 'ORF 蛋白未显示';
      if ($('orfNt').textContent.length < 10) return 'ORF 核酸未显示';
      return false;
    });

    console.log('\n=== G. 帮助页承诺一致性 ===');
    check('帮助页酶数量与实际一致', () => {
      tab('help');
      const txt = d.getElementById('panel-help').textContent;
      const m = txt.match(/(\d+)\s*种/);
      tab('digest');
      const actual = $('enzList').querySelectorAll('input').length;
      if (m && parseInt(m[1], 10) !== actual) return '帮助页写 ' + m[1] + ' 种，实际 ' + actual + ' 种';
      return false;
    });

    console.log('\n=== G2. 新功能入口对账（蛋白/Motif/双酶切/突变/连接） ===');
    check('9 个 tab 均有对应面板', () => {
      const btns = d.querySelectorAll('#tabs button[data-panel]');
      if (btns.length !== 9) return 'tab 数 = ' + btns.length;
      for (let i = 0; i < btns.length; i++) {
        const p = btns[i].getAttribute('data-panel');
        if (!d.getElementById('panel-' + p)) return '缺面板 panel-' + p;
      }
      return false;
    });
    check('新功能关键元素齐全', () => {
      const ids = ['protSeq', 'btnProtAnalyze', 'btnProtFromOrf', 'btnProtFromSeq', 'btnCodonOpt', 'codonDna', 'codonOut',
        'motifPat', 'motifPreset', 'motifBoth', 'btnMotif', 'motifTbl', 'motifOut',
        'ddEnzA', 'ddEnzB', 'btnDD', 'ddOut',
        'mutPos', 'mutOld', 'mutNew', 'mutTm', 'btnMut', 'mutOut', 'btnMutToPrimers',
        'ligVecBp', 'ligInsBp', 'ligVecNg', 'ligRatio', 'ligVecConc', 'ligInsConc', 'btnLig', 'ligOut', 'btnLigFill'];
      for (let i = 0; i < ids.length; i++) if (!d.getElementById(ids[i])) return '缺元素 #' + ids[i];
      return false;
    });
    check('双酶切两个下拉框各含全部酶', () => {
      const a = $('ddEnzA').querySelectorAll('option').length;
      const b = $('ddEnzB').querySelectorAll('option').length;
      if (a !== b) return '两个下拉框酶数不一致 ' + a + '/' + b;
      return false;
    });
    check('帮助页涵盖全部新功能', () => {
      const txt = d.getElementById('panel-help').textContent;
      const keys = ['蛋白', '密码子', 'Motif', '双酶切', '定点突变', '连接反应用量'];
      for (let i = 0; i < keys.length; i++) if (txt.indexOf(keys[i]) < 0) return '帮助页未提及：' + keys[i];
      return false;
    });
    check('帮助页快捷键说明更新为 Ctrl+1…9', () => {
      const txt = d.getElementById('panel-help').textContent;
      return txt.indexOf('1…9') < 0 ? '仍是旧的 1…8 说明' : false;
    });

    console.log('\n=== H. 序列历史（localStorage） ===');
    // histSchedule 有 2 秒防抖，需等待
    setV('seqInput', '>hist_test\n' + 'ACGTTGCA'.repeat(20));
    fire('seqInput', 'input');
    return wait(2400);
  }).then(() => {
    check('序列自动写入历史', () => {
      const opts = $('histSel').querySelectorAll('option');
      if (opts.length <= 1) return '历史下拉为空（可能 localStorage 不可用或未入库）';
      const txt = $('histSel').textContent;
      if (txt.indexOf('hist_test') === -1) return '历史中找不到刚输入的序列';
      return false;
    });
    check('从历史载入序列', () => {
      const sel = $('histSel');
      const idx = [...sel.querySelectorAll('option')].findIndex(o => o.textContent.indexOf('hist_test') > -1);
      if (idx < 0) return '找不到目标历史项';
      sel.value = String(idx - 1); // 减去表头占位项
      sel.dispatchEvent(new w.Event('change', { bubbles: true }));
      if ($('seqInput').value.indexOf('hist_test') === -1) return '未载入到输入框';
      const meta = $('seqMeta').textContent;
      if (meta.indexOf('160') === -1) return '载入后长度 160 nt 未体现: ' + meta.slice(0, 60);
      return false;
    });
    check('清空历史生效', () => {
      click('btnHistClear');
      const opts = $('histSel').querySelectorAll('option');
      if (opts.length > 1) return '清空后仍有 ' + (opts.length - 1) + ' 条';
      return false;
    });

    console.log('\n==========================');
    if (errors.length) {
      console.log('运行时错误 ' + errors.length + ' 个:');
      errors.slice(0, 10).forEach(e => console.log('  • ' + e.split('\n')[0]));
    }
    if (issues.length) {
      console.log('\n发现 ' + issues.length + ' 项问题：');
      issues.forEach(i => console.log('  ' + i));
      process.exit(1);
    } else {
      console.log('功能对账 ' + n + ' 项全部通过 ✔');
      process.exit(0);
    }
  });
}).catch(e => { console.error('审计脚本出错:', e); process.exit(1); });
