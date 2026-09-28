/* 真实 DOM 冒烟测试：用 jsdom 加载 index.html，模拟用户操作，捕获任何运行时异常。
 * 需要 jsdom：NODE_PATH=<node workspace>/node_modules node tools/dom_smoke.js
 */
const path = require('path');
const { JSDOM, VirtualConsole } = require('jsdom');

const root = path.join(__dirname, '..');
const errors = [];

const vc = new VirtualConsole();
vc.on('jsdomError', e => {
  if (/Not implemented/.test(e.message)) return; // scrollTo 等 jsdom 未实现项
  errors.push('jsdomError: ' + (e.stack || e.message));
});
vc.on('error', (...a) => errors.push('console.error: ' + a.join(' ')));

function fakeCtx() {
  const noop = () => {};
  return new Proxy({}, {
    get(t, k) {
      if (k === 'measureText') return () => ({ width: 10 });
      if (k === 'canvas') return null;
      if (k in t) return t[k];
      return noop;
    },
    set(t, k, v) { t[k] = v; return true; }
  });
}

JSDOM.fromFile(path.join(root, 'index.html'), {
  runScripts: 'dangerously',
  resources: 'usable',
  pretendToBeVisual: true,
  virtualConsole: vc,
  beforeParse(win) {
    win.HTMLCanvasElement.prototype.getContext = function () { return fakeCtx(); };
    win.HTMLCanvasElement.prototype.toDataURL = function () { return 'data:image/png;base64,'; };
    win.scrollTo = () => {};
    win.devicePixelRatio = 1;
    win.onerror = (msg, src, line, col, err) => {
      errors.push('window.onerror: ' + msg + ' @' + line + ':' + col + (err && err.stack ? '\n' + err.stack : ''));
    };
    Object.defineProperty(win.Element.prototype, 'clientWidth', { get() { return 900; }, configurable: true });
    Object.defineProperty(win.Element.prototype, 'clientHeight', { get() { return 300; }, configurable: true });
  }
}).then(dom => {
  const win = dom.window;
  return new Promise(resolve => {
    if (win.document.readyState === 'complete') return resolve(dom);
    win.addEventListener('load', () => resolve(dom));
    setTimeout(() => resolve(dom), 4000);
  });
}).then(dom => {
  const win = dom.window, doc = win.document;
  const $ = id => doc.getElementById(id);
  let checks = 0, failed = 0;

  function step(name, fn) {
    const before = errors.length;
    try { fn(); } catch (e) { errors.push('[' + name + '] 抛出异常: ' + (e.stack || e.message)); }
    checks++;
    const added = errors.length - before;
    if (added) { console.log('  [FAIL] ' + name + ' — ' + added + ' 个错误'); failed++; }
    else console.log('  [OK]   ' + name);
  }
  function click(id) {
    const el = $(id);
    if (!el) { errors.push('按钮不存在: ' + id); return; }
    el.dispatchEvent(new win.Event('click', { bubbles: true }));
  }
  function setVal(id, v) { const el = $(id); if (el) el.value = v; }
  function change(id) { const el = $(id); if (el) el.dispatchEvent(new win.Event('change', { bubbles: true })); }

  console.log('\n--- 启动 ---');
  step('页面初始化（自动载入示例序列）', () => {
    const meta = $('seqMeta').textContent;
    if (!/长度/.test(meta)) throw new Error('序列状态条未渲染: ' + meta);
    if ($('ovStats').children.length === 0) throw new Error('概览统计未渲染');
  });
  step('概览面板：统计/组成/分子量已填充', () => {
    if (!/Da/.test($('mwDsDna').textContent)) throw new Error('双链分子量为空');
    if ($('ovComp').innerHTML.length < 20) throw new Error('组成图为空');
  });

  console.log('\n--- 切换每个 Tab ---');
  ['seqops', 'digest', 'primer', 'calc', 'help', 'overview'].forEach(t => {
    step('切换到 ' + t, () => {
      const btn = doc.querySelector('#tabs button[data-panel="' + t + '"]');
      btn.dispatchEvent(new win.Event('click', { bubbles: true }));
      if (!doc.getElementById('panel-' + t).classList.contains('active')) throw new Error('面板未激活');
    });
  });

  console.log('\n--- 序列操作 ---');
  step('反向互补', () => { click('btnRevComp'); if ($('opOut').textContent.length < 10) throw new Error('无输出'); });
  step('互补 / 反向 / 转录 / 逆转录', () => {
    ['btnComp', 'btnRev', 'btnTranscribe', 'btnBackTranscribe'].forEach(click);
  });
  step('六框翻译（+1 ~ −3）', () => {
    ['0', '1', '2', 'r0', 'r1', 'r2'].forEach(f => {
      setVal('frameSel', f); change('frameSel'); click('btnTranslate');
      if ($('protOut').textContent.length < 10) throw new Error('框 ' + f + ' 翻译为空');
    });
  });
  step('ORF 查找并渲染表格', () => {
    click('btnFindOrf');
    const rows = $('orfTbl').querySelectorAll('tbody tr').length;
    if (rows === 0) throw new Error('ORF 表为空');
    if ($('orfProt').textContent.length < 10) throw new Error('ORF 蛋白未显示');
  });
  step('点击 ORF 行的「查看」', () => {
    const b = $('orfTbl').querySelector('button[data-orf]');
    if (b) b.dispatchEvent(new win.Event('click', { bubbles: true }));
  });

  console.log('\n--- 酶切 ---');
  step('载入示例质粒并切到酶切面板', () => {
    doc.querySelector('#tabs button[data-panel="digest"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    setVal('exampleSel', '0'); change('exampleSel');
    if (!$('chkCircular').checked) throw new Error('示例质粒应为环状');
  });
  step('渲染酶列表（51 种）', () => {
    const n = $('enzList').querySelectorAll('input').length;
    if (n < 40) throw new Error('酶列表只有 ' + n + ' 项');
  });
  step('快捷选择：常用 8 种 → 切点与片段表', () => {
    doc.querySelector('[data-quick="common"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    const sites = $('siteTbl').querySelectorAll('tbody tr').length;
    const frags = $('fragTbl').querySelectorAll('tbody tr').length;
    if (sites === 0) throw new Error('切点表为空');
    if (frags === 0) throw new Error('片段表为空');
  });
  step('快捷选择：单切点酶', () => { doc.querySelector('[data-quick="single"]').dispatchEvent(new win.Event('click', { bubbles: true })); });
  step('快捷选择：零切点酶', () => { doc.querySelector('[data-quick="zero"]').dispatchEvent(new win.Event('click', { bubbles: true })); });
  step('快捷选择：pUC19 MCS', () => { doc.querySelector('[data-quick="mcs"]').dispatchEvent(new win.Event('click', { bubbles: true })); });
  step('切换凝胶浓度 / 泳道模式 / 标注', () => {
    ['0.7', '1', '1.5', '2'].forEach(p => {
      setVal('gelPct', p); change('gelPct');
      if ($('gelCanvas').width === 0) throw new Error('凝胶画布宽度为 0');
    });
    $('gelPerEnzyme').checked = false; change('gelPerEnzyme');
    $('gelLabel').checked = false; change('gelLabel');
    $('gelPerEnzyme').checked = true; change('gelPerEnzyme');
  });
  step('酶搜索过滤', () => { setVal('enzSearch', 'eco'); $('enzSearch').dispatchEvent(new win.Event('input', { bubbles: true })); });
  step('清空酶选择', () => { doc.querySelector('[data-quick="none"]').dispatchEvent(new win.Event('click', { bubbles: true })); });
  step('逐个勾选前 5 种酶', () => {
    setVal('enzSearch', ''); $('enzSearch').dispatchEvent(new win.Event('input', { bubbles: true }));
    const inputs = [...$('enzList').querySelectorAll('input')].slice(0, 5);
    inputs.forEach(i => { i.checked = true; i.dispatchEvent(new win.Event('change', { bubbles: true })); });
  });

  console.log('\n--- 引物 & PCR ---');
  step('切到引物面板并自动挑引物', () => {
    doc.querySelector('#tabs button[data-panel="primer"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    click('btnPickFromTpl');
    if (!$('prFwd').value) throw new Error('未挑出正向引物');
    if (!$('prRev').value) throw new Error('未挑出反向引物');
  });
  step('引物参数表渲染', () => {
    const rows = $('prTbl').querySelectorAll('tbody tr').length;
    if (rows !== 2) throw new Error('引物表行数 = ' + rows);
    const tds = $('prTbl').querySelectorAll('tbody tr')[0].querySelectorAll('td');
    if (!/\d/.test(tds[3].textContent)) throw new Error('Tm 未计算');
  });
  step('风险区块渲染', () => { if ($('prRisk').innerHTML.length < 20) throw new Error('风险区为空'); });
  step('PCR 产物预测', () => {
    if (!/bp/.test($('pcrOut').textContent)) throw new Error('产物未输出: ' + $('pcrOut').textContent);
    if ($('pcrSeq').textContent.length < 20) throw new Error('产物序列为空');
  });
  step('PCR 程序表渲染', () => {
    const rows = $('pcrProgram').querySelectorAll('tbody tr').length;
    if (rows < 5) throw new Error('程序表行数 = ' + rows);
  });
  step('产物载入为当前序列', () => { click('btnPcrToInput'); });
  step('手动输入引物 + 改盐浓度', () => {
    setVal('prFwd', 'GTAAAACGACGGCCAGT'); setVal('prRev', 'CAGGAAACAGCTATGAC');
    setVal('prNa', '100'); setVal('prMg', '3'); setVal('prDntp', '0.4'); setVal('prConc', '200');
    change('prNa');
  });
  step('异常输入：非法字符 / 超短引物', () => {
    setVal('prFwd', 'ZZZ'); setVal('prRev', 'A'); change('prFwd');
    setVal('prFwd', ''); setVal('prRev', ''); change('prFwd');
  });

  step('自动引物设计：设计 → 结果表 → 点「使用」填入', () => {
    doc.querySelector('#tabs button[data-panel="primer"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    setVal('exampleSel', '0'); change('exampleSel');
    click('btnDesign');
    const rows = $('dzTbl').querySelectorAll('tbody tr').length;
    if (rows === 0) throw new Error('未设计出候选引物');
    if (/没有满足条件/.test($('dzTbl').textContent)) throw new Error('设计失败: ' + $('dzTbl').textContent.slice(0, 60));
    const useBtn = $('dzTbl').querySelector('button[data-dz]');
    if (!useBtn) throw new Error('缺少「使用」按钮');
    useBtn.dispatchEvent(new win.Event('click', { bubbles: true }));
    if (!$('prFwd').value || !$('prRev').value) throw new Error('引物未填入输入框');
    if (!/bp/.test($('pcrOut').textContent)) throw new Error('填入后未产出 PCR 结果');
  });
  step('自动引物设计：极端参数应给出空态提示而非崩溃', () => {
    setVal('dzMinTm', '95'); setVal('dzMaxTm', '99'); click('btnDesign');
    if (!/没有满足条件/.test($('dzTbl').textContent)) throw new Error('未显示空态提示');
    setVal('dzMinTm', '55'); setVal('dzMaxTm', '65'); click('btnDesign');
  });

  console.log('\n--- 序列比对 ---');
  step('比对：载入演示对 → 统计与视图渲染', () => {
    doc.querySelector('#tabs button[data-panel="align"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    click('btnAlnDemo');
    if (!/一致度/.test($('alnStats').textContent)) throw new Error('未输出一致度: ' + $('alnStats').textContent.slice(0, 50));
    if (!/gap/.test($('alnStats').textContent)) throw new Error('未输出 gap 统计');
    const txt = $('alnOut').textContent;
    if (txt.indexOf('|') === -1) throw new Error('比对视图缺少一致标记');
    if (txt.indexOf('A ') === -1 || txt.indexOf('B ') === -1) throw new Error('比对视图缺少 A/B 行');
  });
  step('比对：当前序列填入与交换', () => {
    setVal('exampleSel', '0'); change('exampleSel');
    click('btnAlnFillA');
    if (!$('alnA').value) throw new Error('未填入 A');
    const before = $('alnA').value;
    click('btnAlnSwap');
    if ($('alnB').value !== before) throw new Error('交换未生效');
  });
  step('比对：超长序列给出提示而非卡死', () => {
    setVal('alnA', 'ACGT'.repeat(800));
    setVal('alnB', 'ACGT'.repeat(800));
    click('btnAlign');
    if (!/过长/.test($('alnStats').textContent)) throw new Error('未提示序列过长');
  });
  step('比对：空输入给出提示', () => {
    setVal('alnA', ''); setVal('alnB', '');
    click('btnAlign');
    if (!/请填写/.test($('alnStats').textContent)) throw new Error('未提示填写序列');
  });
  step('环形图谱：环状显示 / 线性隐藏', () => {
    doc.querySelector('#tabs button[data-panel="digest"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    setVal('exampleSel', '0'); change('exampleSel');   // 环状质粒
    if ($('plasmidCard').style.display === 'none') throw new Error('环状时未显示环形图谱');
    setVal('exampleSel', '1'); change('exampleSel');   // 线性基因组
    if ($('plasmidCard').style.display !== 'none') throw new Error('线性时未隐藏环形图谱');
    setVal('exampleSel', '0'); change('exampleSel');
  });

  console.log('\n--- 虚拟克隆 ---');
  step('克隆：酶切 → 选片段 → 连接 → 重组序列', () => {
    doc.querySelector('#tabs button[data-panel="clone"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    const H = 'GAATTC';
    const vec = 'A'.repeat(500) + H + 'C'.repeat(400) + H + 'G'.repeat(300);
    const ins = H + 'T'.repeat(500) + H;
    setVal('clBackSeq', vec); setVal('clInsSeq', ins);
    $('clBackSeq').dispatchEvent(new win.Event('input', { bubbles: true }));
    $('clInsSeq').dispatchEvent(new win.Event('input', { bubbles: true }));

    const bOpts = $('clBackFrag').querySelectorAll('option').length;
    const iOpts = $('clInsFrag').querySelectorAll('option').length;
    if (bOpts < 2) throw new Error('载体片段下拉未填充（' + bOpts + ' 项）');
    if (iOpts < 2) throw new Error('插入片段下拉未填充（' + iOpts + ' 项）');

    click('btnClone');
    const out = $('clOut').textContent;
    if (out.indexOf('可连接') === -1) throw new Error('未判定为可连接: ' + out.slice(0, 60));
    if (out.indexOf('1312') === -1) throw new Error('重组大小不是 1312 bp: ' + out.slice(0, 80));
    const seqTxt = $('clSeq').textContent.replace(/\n/g, '');
    if (seqTxt.length !== 1312) throw new Error('重组序列长 ' + seqTxt.length + '，期望 1312');
  });
  step('克隆：不兼容组合应给出提示而非静默成功', () => {
    setVal('clInsEnz', 'BamHI');
    $('clInsEnz').dispatchEvent(new win.Event('change', { bubbles: true }));
    click('btnClone');
    const out = $('clOut').textContent;
    if (out.indexOf('不可连接') === -1 && out.indexOf('不匹配') === -1) {
      throw new Error('未提示不兼容: ' + out.slice(0, 80));
    }
    setVal('clInsEnz', 'EcoRI');
    $('clInsEnz').dispatchEvent(new win.Event('change', { bubbles: true }));
  });
  step('克隆：重组序列可载入为当前序列', () => {
    click('btnClone');
    click('btnClToInput');
    if (!$('seqInput').value) throw new Error('未载入到输入框');
    if (!$('chkCircular').checked) throw new Error('重组质粒应自动勾选环状');
    if (!/长度/.test($('seqMeta').textContent)) throw new Error('状态条未更新');
  });

  console.log('\n--- 计算器 ---');
  step('切到计算器面板', () => {
    doc.querySelector('#tabs button[data-panel="calc"]').dispatchEvent(new win.Event('click', { bubbles: true }));
  });
  step('稀释：留空 V1', () => {
    setVal('dC1', '100'); setVal('dV1', ''); setVal('dC2', '10'); setVal('dV2', '100');
    click('btnDilute');
    if (parseFloat($('dV1').value) !== 10) throw new Error('V1 应 = 10，实为 ' + $('dV1').value);
  });
  step('稀释：留空 V2', () => { setVal('dV1', '5'); setVal('dV2', ''); click('btnDilute'); if (!$('dV2').value) throw new Error('V2 未算'); });
  step('稀释：多项留空应给出提示而非崩溃', () => { setVal('dV1', ''); setVal('dV2', ''); click('btnDilute'); });
  step('质量 → mmol', () => { setVal('mMW', '342.3'); setVal('mMass', '34.23'); click('btnMolFromMass'); if (!$('mMol').value) throw new Error('未算出 mmol'); });
  step('mmol → 质量', () => { click('btnMassFromMol'); if (!$('mMass').value) throw new Error('未算出质量'); });
  step('% w/v 换算', () => { setVal('pPct', '1'); setVal('pVol', '100'); setVal('pMW', '342.3'); click('btnPct'); if (!/mg/.test($('molOut').textContent)) throw new Error('未输出'); });
  step('PCR 体系加样表', () => {
    click('btnRx');
    const rows = $('rxTbl').querySelectorAll('tbody tr').length;
    if (rows < 6) throw new Error('加样表行数 = ' + rows);
    if (!/Master Mix/.test($('rxOut').textContent) && !/超/.test($('rxOut').textContent)) throw new Error('无汇总提示');
  });
  step('PCR 体系体积超限时给出警告', () => {
    setVal('rxTplNg', '9999'); click('btnRx');
    if (!/超过|超/.test($('rxOut').textContent)) throw new Error('未提示体积超限');
    setVal('rxTplNg', '50'); click('btnRx');
  });
  step('rpm → ×g', () => { setVal('cfR', '10'); setVal('cfRpm', '12000'); click('btnRpm2g'); if (Math.abs(parseFloat($('cfG').value) - 16099) > 5) throw new Error('×g = ' + $('cfG').value + '，应为 ≈16099'); });
  step('×g → rpm', () => { setVal('cfG', '16000'); click('btnG2rpm'); if (!$('cfRpm').value) throw new Error('未算出 rpm'); });
  step('OD600 菌浓', () => { setVal('od6', '0.6'); setVal('od6dil', '10'); click('btnOd6'); if (!/cells\/mL/.test($('od6Out').textContent)) throw new Error('未输出'); });

  console.log('\n--- 已修 Bug 的回归 ---');
  step('[回归] 翻译结果 HTML 标签不被换行切断', () => {
    doc.querySelector('#tabs button[data-panel="seqops"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    setVal('exampleSel', '0'); change('exampleSel');
    // 找一个足够长的 ORF，保证会跨多行
    setVal('orfMinAa', '60'); click('btnFindOrf');
    const html = $('orfProt').innerHTML;
    const open = (html.match(/<span/g) || []).length;
    const close = (html.match(/<\/span>/g) || []).length;
    if (open === 0) throw new Error('未生成着色 span');
    if (open !== close) throw new Error('span 开闭不平衡: ' + open + ' vs ' + close);
    // 任何一行都不应残留未闭合的 '<'
    const bad = html.split('\n').filter(L => /<[^>]*$/.test(L));
    if (bad.length) throw new Error(bad.length + ' 行残留未闭合标签，例如 ' + JSON.stringify(bad[0].slice(-20)));
  });
  step('[回归] 搜索过滤状态下勾选酶不丢失其它已选酶', () => {
    setVal('exampleSel', '0'); change('exampleSel');
    doc.querySelector('#tabs button[data-panel="digest"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    const before = [...$('enzList').querySelectorAll('input:checked')].map(i => i.value).sort();
    if (before.length < 2) throw new Error('前置条件：应至少选中 2 种酶');
    setVal('enzSearch', 'eco');
    $('enzSearch').dispatchEvent(new win.Event('input', { bubbles: true }));
    const visible = [...$('enzList').querySelectorAll('input')].map(i => i.value);
    // 在过滤结果里切换第一个酶的勾选状态
    const first = $('enzList').querySelector('input');
    first.checked = !first.checked;
    first.dispatchEvent(new win.Event('change', { bubbles: true }));
    setVal('enzSearch', '');
    $('enzSearch').dispatchEvent(new win.Event('input', { bubbles: true }));
    const after = [...$('enzList').querySelectorAll('input:checked')].map(i => i.value);
    // Bug 本质：不在过滤结果里的已选酶不得被丢弃
    const missing = before.filter(n => visible.indexOf(n) === -1 && after.indexOf(n) === -1);
    if (missing.length) throw new Error('丢失了不在过滤结果中的已选酶: ' + missing.join(','));
    // 同时确认被点击的那个酶状态确实翻转了
    const flipped = before.indexOf(first.value) === -1
      ? after.indexOf(first.value) !== -1
      : after.indexOf(first.value) === -1;
    if (!flipped) throw new Error('被点击的酶 ' + first.value + ' 状态未更新');
  });
  step('[回归] 环状序列能识别跨接缝酶切位点', () => {
    setVal('seqInput', '>wrap_test\n' + 'TC' + 'A'.repeat(60) + 'GAAT');
    $('seqInput').dispatchEvent(new win.Event('input', { bubbles: true }));
    $('chkCircular').checked = true; change('chkCircular');
    doc.querySelector('#tabs button[data-panel="digest"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    // 只选 EcoRI
    setVal('enzSearch', 'EcoRI'); $('enzSearch').dispatchEvent(new win.Event('input', { bubbles: true }));
    const inp = $('enzList').querySelector('input');
    if (inp && !inp.checked) { inp.checked = true; inp.dispatchEvent(new win.Event('change', { bubbles: true })); }
    const txt = $('siteTbl').textContent;
    if (txt.indexOf('跨接缝') === -1) throw new Error('切点表未标记跨接缝位点');
  });

  console.log('\n--- 边界情况 ---');
  step('清空序列后各面板不崩溃', () => {
    click('btnClear');
    ['overview', 'seqops', 'digest', 'primer'].forEach(t => {
      doc.querySelector('#tabs button[data-panel="' + t + '"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    });
  });
  step('输入蛋白质序列（应提示而非崩溃）', () => {
    setVal('seqInput', '>prot\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE');
    $('seqInput').dispatchEvent(new win.Event('input', { bubbles: true }));
  });
  // 序列输入有 220ms 防抖，类型判定必须等防抖结束，故放到最后收尾处检查
  step('输入乱码（应提示而非崩溃）', () => {
    setVal('seqInput', '!!!@@@###$$$');
    $('seqInput').dispatchEvent(new win.Event('input', { bubbles: true }));
  });
  step('随机生成 2000 bp', () => { click('btnDemoRandom'); });
  step('GC 窗口参数改为极端值', () => {
    doc.querySelector('#tabs button[data-panel="overview"]').dispatchEvent(new win.Event('click', { bubbles: true }));
    setVal('gcWin', '5'); setVal('gcStep', '1'); change('gcWin');
    setVal('gcWin', '99999'); change('gcWin');
    setVal('gcWin', '100'); setVal('gcStep', '10'); change('gcWin');
  });
  step('浓度互算按钮', () => { click('btnCcFromNg'); click('btnCcFromNm'); });

  // input 事件有 220ms 防抖：类型判定必须等防抖走完，故放到最后收尾
  function wait(ms) { return new Promise(r => setTimeout(r, ms)); }

  return wait(600).then(() => {
    console.log('\n--- 收尾：序列类型判定（需等防抖） ---');
    setVal('seqInput', '>prot2\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE');
    $('seqInput').dispatchEvent(new win.Event('input', { bubbles: true }));
    return wait(400);
  }).then(() => {
    step('[回归] 蛋白序列下序列操作按钮被禁用且有提示', () => {
      doc.querySelector('#tabs button[data-panel="seqops"]').dispatchEvent(new win.Event('click', { bubbles: true }));
      ['btnRevComp', 'btnTranslate', 'btnFindOrf'].forEach(id => {
        if (!$(id).disabled) throw new Error(id + ' 应被禁用');
      });
      if ($('opOut').textContent.indexOf('需要 DNA 或 RNA') === -1) {
        throw new Error('缺少类型提示: ' + $('opOut').textContent.slice(0, 40));
      }
    });
    setVal('seqInput', '>dna2\n' + 'ACGTACGTAC'.repeat(12));
    $('seqInput').dispatchEvent(new win.Event('input', { bubbles: true }));
    return wait(400);
  }).then(() => {
    step('[回归] 切回核酸序列后按钮恢复可用', () => {
      doc.querySelector('#tabs button[data-panel="seqops"]').dispatchEvent(new win.Event('click', { bubbles: true }));
      if ($('btnRevComp').disabled) throw new Error('按钮未恢复可用');
    });
    console.log('\n==========================');
    if (errors.length) {
      console.log('捕获到 ' + errors.length + ' 个运行时错误：\n');
      errors.slice(0, 30).forEach(e => console.log('  • ' + e.split('\n').slice(0, 4).join('\n    ')));
      console.log('\n' + failed + ' / ' + checks + ' 项检查失败 ✘');
      process.exit(1);
    } else {
      console.log(checks + ' 项 DOM 操作检查全部通过，无任何运行时错误 ✔');
      process.exit(0);
    }
  });
}).catch(e => {
  console.error('测试脚本本身出错:', e);
  process.exit(1);
});
