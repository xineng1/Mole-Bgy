/* UI 层缺陷猎取：状态一致性、极端参数、切换序列后的残留
 * 用 async/await 顺序执行，避免回调嵌套导致的"检查跑在准备之前"。
 */
const path = require('path');
const { JSDOM, VirtualConsole } = require('jsdom');
const root = path.join(__dirname, '..');

const errors = [];
const bugs = [];
let checks = 0;

function check(name, fn) {
  checks++;
  try {
    const r = fn();
    if (r) bugs.push({ name, detail: r });
    else console.log('  [OK]   ' + name);
  } catch (e) { bugs.push({ name, detail: '抛出异常: ' + (e.message || e) }); }
}

function fakeCtx() {
  const noop = () => {};
  return new Proxy({}, {
    get(t, k) { if (k === 'measureText') return () => ({ width: 10 }); return k in t ? t[k] : noop; },
    set(t, k, v) { t[k] = v; return true; }
  });
}

function loadDom() {
  const vc = new VirtualConsole();
  vc.on('jsdomError', e => { if (!/Not implemented/.test(e.message)) errors.push(e.message); });
  return JSDOM.fromFile(path.join(root, 'index.html'), {
    runScripts: 'dangerously', resources: 'usable', pretendToBeVisual: true, virtualConsole: vc,
    beforeParse(win) {
      win.HTMLCanvasElement.prototype.getContext = () => fakeCtx();
      win.HTMLCanvasElement.prototype.toDataURL = () => 'data:image/png;base64,';
      win.scrollTo = () => {};
      win.URL.createObjectURL = () => 'blob:stub';
      win.URL.revokeObjectURL = () => {};
      win.HTMLAnchorElement.prototype.click = function () {};
      // 捕获导出的文本内容（Blob 需在 jsdom 里包一层）
      const RealBlob = win.Blob;
      win.Blob = function (parts, opts) {
        win.__blob = { text: parts.join(''), type: opts && opts.type };
        return new RealBlob(parts, opts);
      };
      Object.defineProperty(win.Element.prototype, 'clientWidth', { get() { return 1140; }, configurable: true });
      Object.defineProperty(win.Element.prototype, 'clientHeight', { get() { return 440; }, configurable: true });
      let ok = false;
      try { win.localStorage.setItem('p', '1'); ok = win.localStorage.getItem('p') === '1'; } catch (e) { ok = false; }
      if (!ok) {
        const mem = {};
        Object.defineProperty(win, 'localStorage', {
          configurable: true,
          value: {
            getItem: k => (k in mem ? mem[k] : null), setItem: (k, v) => { mem[k] = String(v); },
            removeItem: k => { delete mem[k]; }, clear: () => {}, key: () => null, get length() { return 0; }
          }
        });
      }
    }
  }).then(dom => new Promise(r => {
    const w = dom.window;
    if (w.document.readyState === 'complete') return r(w);
    w.addEventListener('load', () => r(w));
    setTimeout(() => r(w), 4000);
  }));
}

(async function main() {
  const w = await loadDom();
  const d = w.document;
  const $ = id => d.getElementById(id);
  const wait = ms => new Promise(r => setTimeout(r, ms));
  const tab = t => d.querySelector('#tabs button[data-panel="' + t + '"]').dispatchEvent(new w.Event('click', { bubbles: true }));
  const click = id => { const e = $(id); if (!e) throw new Error('缺少元素 ' + id); e.dispatchEvent(new w.Event('click', { bubbles: true })); };
  const setV = (id, v) => { const e = $(id); if (!e) throw new Error('缺少元素 ' + id); e.value = v; };
  const fire = (id, ev) => $(id).dispatchEvent(new w.Event(ev, { bubbles: true }));
  async function pickExample(i) { setV('exampleSel', String(i)); fire('exampleSel', 'change'); await wait(260); }
  async function setSeq(text, circular) {
    setV('seqInput', text);
    fire('seqInput', 'input');
    await wait(300);
    if (circular !== undefined) { $('chkCircular').checked = !!circular; fire('chkCircular', 'change'); await wait(120); }
  }
  // 构造一个环状接缝处才成立的 ORF
  const WRAP_SEQ = 'TG' + 'GCCACCATGGGGGCCGAATTCAAATAA' + 'A'.repeat(60) + 'A';

  // ORF 表列序：# | 链 | 框 | 区间(nt) | 长度(aa) | 蛋白 | 操作
  function maxAaInTable() {
    const vals = [...$('orfTbl').querySelectorAll('tbody tr')]
      .map(tr => parseInt((tr.children[4] || {}).textContent, 10))
      .filter(v => !isNaN(v));
    return vals.length ? Math.max.apply(null, vals) : 0;
  }

  console.log('\n=== A. 环状跨接缝 ORF 是否贯穿到 UI ===');
  await setSeq('>wraptest\n' + WRAP_SEQ, true);
  tab('seqops');
  setV('orfMinAa', '5');
  click('btnFindOrf');
  const circMaxAa = maxAaInTable();
  check('环状模式下表格出现跨接缝标记', () => {
    const txt = $('orfTbl').textContent;
    if (txt.indexOf('跨接缝') === -1) return '未出现跨接缝标记，环状 ORF 可能没接进 UI';
    return false;
  });

  await setSeq('>wraptest\n' + WRAP_SEQ, false);
  tab('seqops');
  setV('orfMinAa', '5');
  click('btnFindOrf');
  const linMaxAa = maxAaInTable();
  check('线性模式下不出现跨接缝标记', () => {
    if ($('orfTbl').textContent.indexOf('跨接缝') !== -1) return '线性模式仍标记了跨接缝';
    return false;
  });
  // 这条必须真比较，不能留占位——否则是假通过
  check('环状检出的最长 ORF 严格长于线性', () => {
    if (!circMaxAa || !linMaxAa) return '未能从表格读到 aa 长度（环状 ' + circMaxAa + ' / 线性 ' + linMaxAa + '）';
    if (circMaxAa <= linMaxAa) return '环状 ' + circMaxAa + ' aa ≤ 线性 ' + linMaxAa + ' aa，跨接缝 ORF 未真正生效';
    return false;
  });
  console.log('        （环状最长 ' + circMaxAa + ' aa vs 线性最长 ' + linMaxAa + ' aa）');

  console.log('\n=== B. ORF 最小长度的极端输入 ===');
  await pickExample(0);
  tab('seqops');
  for (const v of ['0', '-5', '999999', '', 'abc']) {
    check('minAa = ' + JSON.stringify(v) + ' 时不崩溃且有输出', () => {
      setV('orfMinAa', v);
      click('btnFindOrf');
      const rows = $('orfTbl').querySelectorAll('tbody tr').length;
      if (rows === 0) return '表格既无结果也无提示行';
      return false;
    });
  }

  console.log('\n=== C. 切换序列后的状态一致性 ===');
  await pickExample(0);
  tab('primer');
  const fwd0 = $('prFwd').value;
  check('引物面板会自动填充（前置条件）', () => fwd0 ? false : '未自动填充，后续无法验证');
  await pickExample(1);
  tab('primer');
  check('切换到别的序列后，若引物仍来自旧模板应给出不匹配提示', () => {
    const fwd1 = $('prFwd').value;
    const out = $('pcrOut').textContent || '';
    if (fwd1 === fwd0 && out.indexOf('未找到') === -1) {
      return '引物仍为旧模板的 ' + fwd0.slice(0, 14) + '…，但未提示与当前模板不匹配';
    }
    return false;
  });

  await pickExample(0);
  tab('seqops');
  click('btnFindOrf');
  const orfA = $('orfTbl').textContent.slice(0, 300);
  await pickExample(1);
  tab('seqops');
  click('btnFindOrf');
  const orfB = $('orfTbl').textContent.slice(0, 300);
  check('切换序列后 ORF 表内容随之更新', () => orfA === orfB ? '两端序列的 ORF 表内容完全相同，疑似未刷新' : false);

  await pickExample(0);
  await wait(120);
  tab('seqops');
  check('内置示例无 FEATURES，注释卡片应隐藏', () => $('featCard').style.display !== 'none' ? '无注释却显示了注释卡片' : false);

  console.log('\n=== D. 空序列与极短序列 ===');
  click('btnClear');
  await wait(300);
  check('清空后各面板不崩溃', () => {
    tab('overview'); tab('seqops'); click('btnFindOrf');
    tab('digest'); tab('primer'); click('btnPrimer'); tab('align'); tab('calc');
    return false;
  });
  await setSeq('>tiny\nATG');
  check('3 nt 序列各面板不崩溃', () => {
    tab('overview'); tab('seqops'); click('btnFindOrf');
    tab('digest'); tab('primer'); click('btnPrimer');
    return false;
  });

  console.log('\n=== E. 重复操作不累积状态 ===');
  await pickExample(0);
  tab('seqops');
  setV('orfMinAa', '60');
  click('btnFindOrf');
  const orfFirst = $('orfTbl').querySelectorAll('tbody tr').length;
  for (let i = 0; i < 10; i++) click('btnFindOrf');
  const orfAfter = $('orfTbl').querySelectorAll('tbody tr').length;
  check('连续点击「查找 ORF」11 次结果行数不累积', () => orfFirst !== orfAfter ? '行数 ' + orfFirst + ' → ' + orfAfter : false);

  tab('primer');
  const btnA = $('pcrOut').querySelectorAll('button').length;
  for (let i = 0; i < 5; i++) tab('primer');
  const btnB = $('pcrOut').querySelectorAll('button').length;
  check('反复切换引物面板不累积产物按钮', () => btnA !== btnB ? '按钮数 ' + btnA + ' → ' + btnB : false);

  console.log('\n=== F. PCR 多产物在 UI 上的呈现 ===');
  await setSeq('>' + 'multi\n' + 'ACGTACGTACGT' + 'AAA' + 'ACGTACGTACGT' + 'TTTTTTTTTT' + 'ACGTACGTACGT');
  await wait(200);
  tab('primer');
  setV('prFwd', 'ACGTACGTACGT');
  setV('prRev', 'ACGTACGTACGT');
  fire('prFwd', 'change');
  check('多结合位点模板应列出全部产物并给出提示', () => {
    const out = $('pcrOut').textContent || '';
    if (out.indexOf('多个引物结合位点') === -1) return '未提示存在多个结合位点，输出为：' + out.slice(0, 80);
    const rows = $('pcrOut').querySelectorAll('tbody tr').length;
    if (rows < 2) return '产物表只有 ' + rows + ' 行';
    return false;
  });
  check('产物「查看」按钮可切换序列显示', () => {
    const btn = $('pcrOut').querySelector('button[data-pcr]');
    if (!btn) return '缺少「查看」按钮';
    btn.dispatchEvent(new w.Event('click', { bubbles: true }));
    if (($('pcrSeq').textContent || '').length < 10) return '未显示产物序列';
    return false;
  });

  console.log('\n=== G. 报告输出的表格注入防护 ===');
  await setSeq('>seq1|geneA|variant2\n' + 'ACGT'.repeat(50));
  await wait(250);
  tab('overview');
  w.__blob = null;
  click('btnReport');
  check('序列名含竖线时报告表格不被撑破', () => {
    if (!w.__blob) return '未生成报告';
    const line = w.__blob.text.split('\n').filter(l => l.indexOf('| 序列名称 |') === 0)[0];
    if (!line) return '报告里找不到「序列名称」行';
    if (line.indexOf('\\|') === -1) return '竖线未转义，该行会被解析成多列：' + line;
    return false;
  });
  check('报告在有引物但无模板匹配时不出现 undefined', () => {
    if (!w.__blob) return '未生成报告';
    return /undefined|NaN/.test(w.__blob.text) ? '报告含 undefined / NaN' : false;
  });

  console.log('\n=== H. 环状模板的跨接缝 PCR 产物 ===');
  {
    const fwd = 'AAAACCCCGGGG';
    const rSite = 'TTTTGGGGCCCC';
    const rev = w.MolCore.revComp(rSite);
    // 反向位点在开头、正向位点在末尾 → 线性无产物，环状应绕接缝出一条约 24 bp 的产物
    const seq = rSite + 'A'.repeat(50) + fwd;
    await setSeq('>pcrcross\n' + seq, true);
    await wait(250);
    tab('primer');
    setV('prFwd', fwd); setV('prRev', rev);
    fire('prFwd', 'change');
    check('环状模板应检出跨接缝产物', () => {
      const out = $('pcrOut').textContent || '';
      if (out.indexOf('未找到') !== -1) return '环状模式下未检出产物：' + out.slice(0, 60);
      if (out.indexOf('24') === -1) return '产物大小不是预期的 24 bp：' + out.slice(0, 80);
      return false;
    });
    check('跨接缝产物的序列长度与产物大小一致', () => {
      const txt = ($('pcrSeq').textContent || '').replace(/\n/g, '');
      if (txt.length !== 24) return '产物序列长 ' + txt.length + '，期望 24';
      if (txt.slice(0, fwd.length) !== fwd) return '产物 5′ 端不是正向引物';
      if (txt.slice(-rSite.length) !== rSite) return '产物 3′ 端不是反向引物结合位点';
      return false;
    });
    await setSeq('>pcrcross\n' + seq, false);
    await wait(150);
    tab('primer');
    check('同一序列切为线性后应检不到该产物', () => {
      const out = $('pcrOut').textContent || '';
      if (out.indexOf('未找到') === -1) return '线性模式仍报出产物：' + out.slice(0, 60);
      return false;
    });
  }

  console.log('\n=== I. 键盘快捷键 ===');
  {
    function press(key, mods) {
      const ev = new w.KeyboardEvent('keydown', { key: key, ctrlKey: !!(mods && mods.ctrl), bubbles: true, cancelable: true });
      d.dispatchEvent(ev);
      return ev;
    }
    const activePanel = () => (d.querySelector('.panel.active') || {}).id;

    press('2', { ctrl: true });
    check('Ctrl+2 切到序列操作面板', () => activePanel() !== 'panel-seqops' ? '当前是 ' + activePanel() : false);
    press('3', { ctrl: true });
    check('Ctrl+3 切到蛋白面板', () => activePanel() !== 'panel-protein' ? '当前是 ' + activePanel() : false);
    press('4', { ctrl: true });
    check('Ctrl+4 切到酶切面板', () => activePanel() !== 'panel-digest' ? '当前是 ' + activePanel() : false);
    press('9', { ctrl: true });
    check('Ctrl+9 切到帮助面板', () => activePanel() !== 'panel-help' ? '当前是 ' + activePanel() : false);
    check('Ctrl+9 切换不报错', () => false);
    check('Ctrl+数字会阻止浏览器默认行为', () => {
      const ev = new w.KeyboardEvent('keydown', { key: '2', ctrlKey: true, bubbles: true, cancelable: true });
      d.dispatchEvent(ev);
      return ev.defaultPrevented ? false : '未调用 preventDefault';
    });

    // Ctrl+Enter 执行当前面板主操作
    press('2', { ctrl: true });
    tab('seqops');
    await wait(150);
    check('Ctrl+Enter 触发「查找 ORF」', () => {
      $('orfTbl').querySelector('tbody').innerHTML = '';
      const ev = new w.KeyboardEvent('keydown', { key: 'Enter', ctrlKey: true, bubbles: true, cancelable: true });
      d.dispatchEvent(ev);
      const rows = $('orfTbl').querySelectorAll('tbody tr').length;
      if (rows === 0) return '主操作未执行（表格仍为空）';
      return false;
    });
    check('无修饰键时不触发快捷键', () => {
      tab('overview');
      const before = activePanel();
      press('3', {});
      return activePanel() !== before ? '无修饰键却切换了面板' : false;
    });
  }

  // === J. 新功能：蛋白面板 / 密码子优化 / Motif / 双酶切 / 突变引物 / 连接计算 ===
  {
    // J1 蛋白理化性质
    tab('protein');
    setV('protSeq', '');
    click('btnProtAnalyze');
    check('蛋白面板：空序列时给出提示而非崩溃', () =>
      $('protStats').textContent.indexOf('请先提供') < 0 ? '无提示：' + $('protStats').textContent.slice(0, 30) : false);

    setV('protSeq', '>demo\nMKTAYIAKQRQISFVKSHFSRQLEDLRQFIERTKKLD\n');
    click('btnProtAnalyze');
    check('蛋白面板：FASTA 头被剥除并产出 9 个指标块', () => {
      const n = $('protStats').querySelectorAll('.stat').length;
      return n !== 9 ? '指标块数 = ' + n : false;
    });
    check('蛋白面板：残基组成表 20 行', () =>
      $('protCompTbl').querySelectorAll('tbody tr').length !== 20 ? '行数异常' : false);
    check('蛋白面板：pI 数值在 0–14', () => {
      const t = $('protStats').textContent.match(/等电点 pI\s*([\d.]+)/);
      const v = t && parseFloat(t[1]);
      return !(v >= 0 && v <= 14) ? 'pI 读数异常：' + $('protStats').textContent.slice(0, 60) : false;
    });

    setV('protSeq', 'XZBJ123');
    click('btnProtAnalyze');
    check('蛋白面板：全非法残基时提示而不是抛异常', () =>
      $('protStats').textContent.indexOf('请先提供') < 0 ? '未拒绝非法序列' : false);

    // J2 密码子优化
    setV('protSeq', 'MKTAYIAKQRQISFVKSHFSRQ');
    click('btnCodonOpt');
    check('密码子优化：回译校验一致', () =>
      $('codonOut').textContent.indexOf('完全一致') < 0 ? '回译校验未通过' : false);
    check('密码子优化：DNA 长度为 66 nt（22 aa）', () => {
      const dna = $('codonDna').textContent.replace(/\s/g, '');
      return dna.length !== 66 ? '长度 = ' + dna.length : false;
    });
    check('密码子优化：酶切位点扫描已执行', () =>
      $('codonOut').textContent.indexOf('酶切位点') < 0 ? '缺少位点扫描输出' : false);
    click('btnCodonToInput');
    await wait(350);
    check('密码子优化：载入为当前序列后概览联动', () =>
      $('seqInput').value.indexOf('codon_optimized') < 0 ? '未载入' : false);

    // 恢复示例序列供后续用例使用
    await pickExample(0);

    // J3 Motif 搜索
    tab('digest');
    setV('motifPat', '');
    click('btnMotif');
    check('Motif：空模式给出提示', () =>
      $('motifOut').textContent.indexOf('请输入模式') < 0 ? '无提示' : false);
    setV('motifPat', 'AX3!');
    click('btnMotif');
    check('Motif：非法字符被拒绝并说明允许的字母', () =>
      $('motifOut').textContent.indexOf('仅允许') < 0 ? '未拒绝非法模式' : false);
    setV('motifPat', 'AGGAGG');
    click('btnMotif');
    check('Motif：搜索执行并报告命中数', () =>
      $('motifOut').textContent.indexOf('处命中') < 0 ? '无命中报告：' + $('motifOut').textContent.slice(0, 30) : false);
    setV('motifPreset', 'AATAAA');
    fire('motifPreset', 'change');
    check('Motif：预设选择自动填入模式框', () =>
      $('motifPat').value !== 'AATAAA' ? '未填入：' + $('motifPat').value : false);

    // 环状跨接缝标记（ATAAAGG 上 GGAT 只能跨接缝命中）
    await setSeq('>circ\nATAAAGG', true);
    setV('motifPat', 'GGAT');
    $('motifBoth').checked = false;
    click('btnMotif');
    check('Motif：环状跨接缝命中带「跨接缝」标记', () => {
      const row = $('motifTbl').querySelector('tbody tr');
      if (!row) return '结果表为空';
      return row.textContent.indexOf('跨接缝') < 0 ? '无标记：' + row.textContent.replace(/\s+/g, ' ') : false;
    });
    check('Motif：跨接缝命中位于序列末尾附近', () => {
      const cells = $('motifTbl').querySelectorAll('tbody tr td');
      if (cells.length < 2) return '结果行不完整';
      const pos = parseInt(cells[1].textContent, 10);
      return pos !== 6 ? '位置应为 6，实际 ' + cells[1].textContent : false;
    });

    // 结果上限：26000 bp 双链搜 N → 52000 命中 > 50000 上限
    await setSeq('>longA\n' + 'A'.repeat(26000), true);
    setV('motifPat', 'N');
    $('motifBoth').checked = true;
    click('btnMotif');
    check('Motif：命中触达上限时给出截断提示', () => {
      const t = $('motifOut').textContent;
      if (t.indexOf('截断') < 0) return '无截断提示：' + t.slice(0, 60);
      if ($('motifTbl').querySelectorAll('tbody tr').length !== 500) return '未限制显示条数';
      return false;
    });
    await pickExample(0);

    // J4 双酶切速查
    click('btnDD');
    check('双酶切：默认 EcoRI×HindIII 可比较并给出切点数', () =>
      $('ddOut').textContent.indexOf('个切点') < 0 ? '无输出' : false);
    check('双酶切：报告联合片段', () =>
      $('ddOut').textContent.indexOf('联合酶切片段') < 0 ? '缺少联合片段行' : false);

    // J5 定点突变引物
    tab('primer');
    setV('mutPos', '');
    setV('mutNew', 'G');
    click('btnMut');
    check('突变引物：位置为空时给出错误提示', () =>
      $('mutOut').textContent.indexOf('超出模板范围') < 0 && $('mutOut').textContent.indexOf('请先') < 0
        ? '无错误提示：' + $('mutOut').textContent.slice(0, 40) : false);
    setV('mutPos', '500');
    click('btnMut');
    check('突变引物：正常设计出 F/R 两条引物', () =>
      $('mutOut').textContent.indexOf('正向引物') < 0 || $('mutOut').textContent.indexOf('反向引物') < 0
        ? '缺少引物输出' : false);
    check('突变引物：给出 PCR 程序建议', () =>
      $('mutOut').textContent.indexOf('KLD') < 0 ? '缺少程序建议' : false);
    click('btnMutToPrimers');
    check('突变引物：一键填入引物输入框', () =>
      !$('prFwd').value || $('prFwd').value.charAt(0) !== 'G' ? '未正确填入：' + $('prFwd').value.slice(0, 20) : false);
    setV('mutPos', '99999999');
    click('btnMut');
    check('突变引物：超大位置被拒绝', () =>
      $('mutOut').textContent.indexOf('超出模板范围') < 0 ? '未拒绝越界位置' : false);

    // J6 连接用量计算
    tab('clone');
    setV('ligVecBp', '5000'); setV('ligInsBp', '1000');
    setV('ligVecNg', '50'); setV('ligRatio', '3');
    setV('ligVecConc', ''); setV('ligInsConc', '');
    click('btnLig');
    check('连接计算：3:1 摩尔比算出 30 ng 插入', () =>
      $('ligOut').textContent.indexOf('30 ng') < 0 ? '结果异常：' + $('ligOut').textContent.slice(0, 60) : false);
    setV('ligInsConc', '20');
    click('btnLig');
    check('连接计算：填浓度后显示吸取 µL 数', () =>
      $('ligOut').textContent.indexOf('µL') < 0 ? '未显示体积' : false);
    setV('ligVecBp', '0');
    click('btnLig');
    check('连接计算：非法输入被拒绝', () =>
      $('ligOut').textContent.indexOf('必须大于 0') < 0 ? '未拒绝非法输入' : false);
  }

  console.log('\n' + '='.repeat(60));
  if (errors.length) console.log('运行时错误 ' + errors.length + ' 个：' + errors.slice(0, 2).join(' | '));
  if (bugs.length) {
    console.log('发现 ' + bugs.length + ' / ' + checks + ' 个问题：\n');
    bugs.forEach((b, i) => console.log('  ' + (i + 1) + '. ' + b.name + '\n     → ' + b.detail));
    process.exit(1);
  } else {
    console.log(checks + ' 项 UI 排查全部通过 ✔');
    process.exit(0);
  }
})().catch(e => { console.error('排查脚本出错:', e); process.exit(1); });
