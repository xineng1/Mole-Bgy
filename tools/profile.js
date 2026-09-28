/* 性能剖析：在 jsdom 中模拟真实操作，定位耗时热点 */
const path = require('path');
const { JSDOM, VirtualConsole } = require('jsdom');
const root = path.join(__dirname, '..');

const vc = new VirtualConsole();
const errors = [];
vc.on('jsdomError', e => { if (!/Not implemented/.test(e.message)) errors.push(e.message); });

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
    win.scrollTo = () => {};
    // jsdom 未实现 download 相关 API，stub 掉以免污染性能基线
    win.URL.createObjectURL = () => 'blob:stub';
    win.URL.revokeObjectURL = () => {};
    win.HTMLAnchorElement.prototype.click = function () {};
    Object.defineProperty(win.Element.prototype, 'clientWidth', { get() { return 1140; }, configurable: true });
    Object.defineProperty(win.Element.prototype, 'clientHeight', { get() { return 440; }, configurable: true });
    let ok = false;
    try { win.localStorage.setItem('p', '1'); ok = win.localStorage.getItem('p') === '1'; } catch (e) { ok = false; }
    if (!ok) {
      const mem = {};
      Object.defineProperty(win, 'localStorage', {
        configurable: true,
        value: { getItem: k => (k in mem ? mem[k] : null), setItem: (k, v) => { mem[k] = String(v); }, removeItem: k => { delete mem[k]; }, clear: () => {}, key: () => null, get length() { return 0; } }
      });
    }
  }
}).then(dom => new Promise(r => {
  const w = dom.window;
  if (w.document.readyState === 'complete') return r(w);
  w.addEventListener('load', () => r(w));
  setTimeout(() => r(w), 4000);
})).then(w => {
  const d = w.document, $ = id => d.getElementById(id);
  // 单次测量在 jsdom 下噪声很大（GC / JIT 预热），对每个操作多轮采样取最小值，
  // 最小值最接近该操作的"纯粹成本"，用于优化前后对比才可信。
  const samples = new Map();
  function t(label, fn) {
    const t0 = Date.now();
    fn();
    const ms = Date.now() - t0;
    const arr = samples.get(label) || [];
    arr.push(ms);
    samples.set(label, arr);
    return ms;
  }
  function report() {
    const rows = [];
    samples.forEach((arr, label) => rows.push({ label, min: Math.min.apply(null, arr), n: arr.length, arr }));
    rows.sort((a, b) => b.min - a.min);
    const fmt = v => String(v).padStart(6);
    console.log('\n  ' + fmt('最小') + '  ' + '样本数  操作（多轮取最小，越小越接近真实成本）');
    console.log('  ' + '-'.repeat(76));
    rows.forEach(r => console.log('  ' + fmt(r.min) + ' ms   ×' + String(r.n).padEnd(5) + '  ' + r.label + '   ' + JSON.stringify(r.arr)));
    const slow = rows.filter(r => r.min >= 200);
    if (slow.length) {
      console.log('\n  仍超过 200 ms 的操作：');
      slow.forEach(r => console.log('    • ' + fmt(r.min) + ' ms  ' + r.label));
    } else {
      console.log('\n  所有操作均已低于 200 ms 阈值 ✔');
    }
  }
  function tab(name) {
    d.querySelector('#tabs button[data-panel="' + name + '"]').dispatchEvent(new w.Event('click', { bubbles: true }));
  }
  function setSeq(text) {
    $('seqInput').value = text;
    $('seqInput').dispatchEvent(new w.Event('input', { bubbles: true }));
  }
  const wait = ms => new Promise(r => setTimeout(r, ms));

  // 参照：核心层单次调用的裸成本
  function mk(n) { return w.MolCore.randomSeq(n, 50); }
  const big = mk(120000);
  const tb = Date.now(); w.MolCore.findORFs(big, 30, true);
  console.log('参照：核心层 findORFs(120 kb) 单次 = ' + (Date.now() - tb) + ' ms（每次切面板都会重算）');

  console.log('\n=== 场景 A：载入 120 kb 序列 ===');
  const tA = Date.now();
  setSeq('>perf120k\n' + big);
  return wait(700).then(() => {
    // 载入大序列代价高，只采一次；其中固定含 220 ms 防抖
    samples.set('输入 120 kb（含 220 ms 防抖）', [Date.now() - tA]);

    // 统计 findORFs 的真实调用次数：缓存生效后，反复切面板不应重复触发六框扫描
    const origFindORFs = w.MolCore.findORFs;
    let findCalls = 0;
    w.MolCore.findORFs = function () { findCalls++; return origFindORFs.apply(this, arguments); };

    console.log('\n=== 场景 B/C/D：多轮采样（每项 3 次取最小） ===');
    for (let round = 0; round < 3; round++) {
      ['seqops', 'digest', 'primer', 'overview'].forEach(n => t('切到 ' + n, () => tab(n)));

      tab('digest');
      t('勾选一个酶后重绘', () => {
        const inp = $('enzList').querySelector('input');
        inp.checked = !inp.checked;
        inp.dispatchEvent(new w.Event('change', { bubbles: true }));
      });
      t('切换凝胶浓度后重绘', () => {
        $('gelPct').value = $('gelPct').value === '1.5' ? '1' : '1.5';
        $('gelPct').dispatchEvent(new w.Event('change', { bubbles: true }));
      });
      t('酶搜索输入（重建列表）', () => {
        $('enzSearch').value = 'ec';
        $('enzSearch').dispatchEvent(new w.Event('input', { bubbles: true }));
      });

      tab('primer');
      t('改一个引物碱基（重算风险区）', () => {
        $('prFwd').value = $('prFwd').value.slice(0, -1) + 'A';
        $('prFwd').dispatchEvent(new w.Event('change', { bubbles: true }));
      });
      t('自动引物设计（默认全长）', () => {
        $('btnDesign').dispatchEvent(new w.Event('click', { bubbles: true }));
      });
    }

    console.log('\n=== 场景 E：比对 1500×1500 ===');
    tab('align');
    const a = w.MolCore.randomSeq(1500), b = w.MolCore.randomSeq(1500);
    t('填入两条 1500 nt 序列比对面板', () => { $('alnA').value = a; $('alnB').value = b; });
    t('执行比对', () => { $('btnAlign').dispatchEvent(new w.Event('click', { bubbles: true })); });

    console.log('\n=== 场景 F：导出报告 ===');
    tab('overview');
    t('生成并导出 Markdown 报告', () => { $('btnReport').dispatchEvent(new w.Event('click', { bubbles: true })); });

    console.log('\n=== 缓存有效性 ===');
    const perCall = w.MolCore.findORFs; // 保持包装
    console.log('  3 轮切换 + 3 次勾选酶 + 3 次改凝胶浓度 + 3 次引物设计，');
    console.log('  在这期间 findORFs 实际被调用 ' + findCalls + ' 次。');
    console.log('  由于序列未变，理想值应为 2 次（minAa=30 与 60 各一次），其余全部命中缓存。');
    console.log('  若未加缓存，每次切到酶切面板要算 2 次、概览 1 次，3 轮约 9 次以上。');

    console.log('\n================ 汇总 ================');
    report();
    if (errors.length) console.log('\n运行时错误: ' + errors.length);
    process.exit(0);
  });
}).catch(e => { console.error('剖析脚本出错:', e); process.exit(1); });
