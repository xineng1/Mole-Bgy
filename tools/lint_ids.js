/* 静态检查：JS 语法 + HTML/JS 之间的 DOM id 引用一致性 */
const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');
const root = path.join(__dirname, '..');

let fail = 0;

// 1. JS 语法检查
['core.js', 'enzymes.js', 'examples.js', 'app.js'].forEach(f => {
  try {
    execFileSync(process.execPath, ['--check', path.join(root, 'assets', f)], { stdio: 'pipe' });
    console.log('  [OK]   语法 ' + f);
  } catch (e) {
    console.log('  [FAIL] 语法 ' + f + '\n' + e.stderr.toString());
    fail++;
  }
});

// 2. id 一致性
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const app = fs.readFileSync(path.join(root, 'assets', 'app.js'), 'utf8');

const htmlIds = new Set();
html.replace(/id="([^"]+)"/g, (m, id) => { htmlIds.add(id); return m; });

const usedIds = new Set();
app.replace(/\$\('([^']+)'\)/g, (m, id) => { usedIds.add(id); return m; });
// querySelector 里用到的 id 选择器
app.replace(/querySelector(?:All)?\('#([A-Za-z0-9_-]+)/g, (m, id) => { usedIds.add(id); return m; });

const missing = [...usedIds].filter(id => !htmlIds.has(id));
if (missing.length) {
  console.log('  [FAIL] app.js 引用了 HTML 中不存在的 id: ' + missing.join(', '));
  fail++;
} else {
  console.log('  [OK]   app.js 引用的 ' + usedIds.size + ' 个 id 全部存在于 index.html');
}

const unused = [...htmlIds].filter(id => !usedIds.has(id) && !/^panel-/.test(id) && id !== 'tabs');
if (unused.length) console.log('  [info] HTML 中未被 app.js 直接引用的 id: ' + unused.join(', '));

// 3. 脚本引用检查
const scripts = [...html.matchAll(/<script src="([^"]+)"/g)].map(m => m[1]);
scripts.forEach(s => {
  const p = path.join(root, s);
  const exists = fs.existsSync(p);
  console.log((exists ? '  [OK]   ' : '  [FAIL] ') + '脚本引用 ' + s);
  if (!exists) fail++;
});
const links = [...html.matchAll(/<link[^>]+href="([^"]+)"/g)].map(m => m[1]);
links.forEach(s => {
  const exists = fs.existsSync(path.join(root, s));
  console.log((exists ? '  [OK]   ' : '  [FAIL] ') + '样式引用 ' + s);
  if (!exists) fail++;
});

// 4. tab 面板完整性
const tabs = [...html.matchAll(/data-panel="([^"]+)"/g)].map(m => m[1]);
tabs.forEach(t => {
  const okp = htmlIds.has('panel-' + t);
  console.log((okp ? '  [OK]   ' : '  [FAIL] ') + 'tab 面板 panel-' + t);
  if (!okp) fail++;
});

// 5. 酶库完整性：每个酶的切点偏移必须落在识别序列范围内
const vm = require('vm');
function load(f) {
  const ctx = { window: {}, console };
  vm.createContext(ctx);
  vm.runInContext(fs.readFileSync(path.join(root, 'assets', f), 'utf8'), ctx, { filename: f });
  return ctx.window;
}
const E = load('enzymes.js').MolEnzymes;
let bad = [];
E.list.forEach(e => {
  if (e.cutOffset < 0 || e.cutOffset > e.siteClean.length) bad.push(e.name);
  if (!/^[ACGTRYSWKMBDHVN]+$/.test(e.siteClean)) bad.push(e.name + '(非法字符)');
});
console.log(bad.length ? '  [FAIL] 酶库异常: ' + bad.join(', ') : '  [OK]   酶库 ' + E.list.length + ' 种，切点偏移与字符全部合法');
if (bad.length) fail++;
const dup = E.names.filter((n, i) => E.names.indexOf(n) !== i);
if (dup.length) { console.log('  [FAIL] 重复酶名: ' + dup.join(',')); fail++; }

console.log('\n' + (fail === 0 ? '静态检查全部通过 ✔' : fail + ' 项失败 ✘'));
process.exit(fail === 0 ? 0 : 1);
