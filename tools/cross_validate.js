/* 算法交叉验证：用独立的朴素实现去对拍核心函数，找出实现层难以自测的错误。
 * 朴素实现只求"直觉上显然正确"，不追求性能，用来当参照系。
 */
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const root = path.join(__dirname, '..');

function load(f) {
  const ctx = { window: {}, console };
  vm.createContext(ctx);
  vm.runInContext(fs.readFileSync(path.join(root, 'assets', f), 'utf8'), ctx, { filename: f });
  return ctx.window;
}
const C = load('core.js').MolCore;
const E = load('enzymes.js').MolEnzymes;
const EX = load('examples.js').MolExamples;

const fails = [];
function eq(name, got, want, extra) {
  const g = JSON.stringify(got), w = JSON.stringify(want);
  if (g !== w) fails.push({ name, got: g, want: w, extra: extra || '' });
}

function randSeq(n, seed) {
  let s = '', x = seed || 12345;
  const B = 'ACGT';
  for (let i = 0; i < n; i++) { x = (x * 1103515245 + 12345) & 0x7fffffff; s += B[x % 4]; }
  return s;
}

/* ---------- 朴素实现 ---------- */
const IUPAC_MAP = {
  A: 'A', C: 'C', G: 'G', T: 'T', U: 'T',
  R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC',
  B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT'
};
function baseMatches(pattern, base) { return (IUPAC_MAP[pattern] || pattern).indexOf(base) !== -1; }

/** 朴素酶切：逐位置试匹配，环状时下标取模 */
function naiveSites(seq, enzyme, circular) {
  const site = enzyme.site.replace('^', '');
  const cutOff = enzyme.site.indexOf('^') === -1 ? Math.floor(site.length / 2) : enzyme.site.indexOf('^');
  const n = seq.length;
  const limit = circular ? n : n - site.length + 1;
  const cuts = [];
  for (let i = 0; i < limit; i++) {
    let ok = true;
    for (let j = 0; j < site.length; j++) {
      if (!baseMatches(site[j], seq[(i + j) % n])) { ok = false; break; }
    }
    if (ok) cuts.push((i + cutOff) % n);
  }
  return cuts.filter(p => p > 0 && p < n).sort((a, b) => a - b);
}

function naiveFragments(len, cuts, circular) {
  const u = [...new Set(cuts)].sort((a, b) => a - b);
  if (!u.length) return [len];
  if (circular) {
    if (u.length === 1) return [len];
    const fr = [];
    for (let i = 1; i < u.length; i++) fr.push(u[i] - u[i - 1]);
    fr.push(len + u[0] - u[u.length - 1]);
    return fr.sort((a, b) => b - a);
  }
  const fr = [u[0]];
  for (let i = 1; i < u.length; i++) fr.push(u[i] - u[i - 1]);
  fr.push(len - u[u.length - 1]);
  return fr.sort((a, b) => b - a);
}

/** 朴素 ORF：逐条链、逐个框，暴力找 ATG…终止子 */
function naiveORFs(seq, minAa) {
  const CODON = C.CODON_TABLE;
  const n = seq.length;
  const out = [];
  [['+', seq], ['-', C.revComp(seq)]].forEach(([strand, s]) => {
    for (let f = 0; f < 3; f++) {
      for (let i = f; i + 3 <= s.length; i += 3) {
        if (s.substr(i, 3) !== 'ATG') continue;
        for (let k = i; k + 3 <= s.length; k += 3) {
          const aa = CODON[s.substr(k, 3)];
          if (aa === '*') {
            const len = (k - i) / 3;
            // core 的 frame 是 1-based，这里对齐
            if (len >= minAa) out.push(strand + (f + 1) + ':' + len);
            break;
          }
        }
      }
    }
  });
  return out.sort();
}

/** 朴素 PCR 产物 */
function naiveProducts(t, fwd, rev, circular) {
  const rSite = C.revComp(rev);
  const findAll = (hay, needle) => {
    const r = []; let i = hay.indexOf(needle);
    while (i !== -1) { r.push(i); i = hay.indexOf(needle, i + 1); }
    return r;
  };
  const fs = findAll(t, fwd), rs = findAll(t, rSite), n = t.length;
  const out = [];
  fs.forEach(a => rs.forEach(b => {
    if (b > a) out.push(b + rSite.length - a);
    else if (circular && b !== a) {
      const size = (n - a) + b + rSite.length;
      if (size <= n) out.push(size);
    }
  }));
  return out.sort((x, y) => x - y);
}

/* ---------- 1. 酶切 ---------- */
console.log('=== 1. 酶切位点与片段（对拍朴素实现）===');
{
  let cases = 0, mismatch = 0, firstBad = null;
  const enzymes = ['EcoRI', 'BamHI', 'HindIII', 'PstI', 'TaqI', 'AluI', 'MboI', 'AvaI', 'HhaI', 'NotI', 'SmaI'];
  for (let seed = 1; seed <= 40; seed++) {
    const len = 60 + (seed * 37) % 400;
    const seq = randSeq(len, seed * 991);
    for (const name of enzymes) {
      const en = E.byName(name);
      for (const circular of [false, true]) {
        cases++;
        const mine = C.findSites(seq, en, circular).map(h => h.cut).filter(p => p > 0 && p < len).sort((a, b) => a - b);
        const want = naiveSites(seq, en, circular);
        if (JSON.stringify(mine) !== JSON.stringify(want)) {
          mismatch++;
          if (!firstBad) firstBad = { seed, name, circular, len, mine: mine.slice(0, 8), want: want.slice(0, 8) };
        }
        // 片段长度也要一致
        const myFrag = C.fragmentsFromCuts(len, mine, circular);
        const wantFrag = naiveFragments(len, want, circular);
        if (JSON.stringify(myFrag) !== JSON.stringify(wantFrag) && !firstBad) {
          mismatch++;
          firstBad = { seed, name, circular, len, mine: myFrag, want: wantFrag, note: '片段不一致' };
        }
      }
    }
  }
  eq('酶切对拍（' + cases + ' 组）', mismatch === 0 ? 'ok' : firstBad, 'ok',
    mismatch ? '' : '');
  console.log('  对比 ' + cases + ' 组，不一致 ' + mismatch + ' 组');
  if (firstBad) console.log('  首个不一致：' + JSON.stringify(firstBad));
}

/* ---------- 2. ORF ---------- */
console.log('\n=== 2. ORF 六框扫描（对拍朴素实现）===');
{
  let cases = 0, mismatch = 0, firstBad = null;
  for (let seed = 1; seed <= 25; seed++) {
    const seq = randSeq(150 + seed * 13, seed * 7717);
    // 混入一些 ATG/终止子，让 ORF 更容易出现
    const injected = seq.slice(0, 50) + 'ATGGCCATTGTAATGGGCCGCTAA' + seq.slice(50);
    for (const minAa of [1, 3, 8]) {
      cases++;
      const mine = C.findORFs(injected, minAa, true).map(o => o.strand + o.frame + ':' + o.aaLength).sort();
      const want = naiveORFs(injected, minAa);
      if (JSON.stringify(mine) !== JSON.stringify(want)) {
        mismatch++;
        if (!firstBad) firstBad = { seed, minAa, mine: mine.slice(0, 6), want: want.slice(0, 6) };
      }
    }
  }
  console.log('  对比 ' + cases + ' 组，不一致 ' + mismatch + ' 组');
  if (firstBad) console.log('  首个不一致：' + JSON.stringify(firstBad));
  eq('ORF 对拍', mismatch, 0);
}

/* ---------- 3. 环状 ORF 不越过一圈 ---------- */
console.log('\n=== 3. 环状 ORF 的规模与去重 ===');
{
  let bad = 0, maxLen = 0, dup = 0;
  for (let seed = 1; seed <= 12; seed++) {
    const n = 200 + seed * 7;
    const seq = randSeq(n, seed * 331);
    const orfs = C.findORFsCircular(seq, 2, true);
    const seen = new Set();
    orfs.forEach(o => {
      const len = o.nt.length;
      if (len > n) bad++;
      maxLen = Math.max(maxLen, len);
      const key = o.strand + ':' + o.ntStart + ':' + o.ntEnd;
      if (seen.has(key)) dup++;
      seen.add(key);
      if (!o.nt || o.nt.length !== (o.aaLength + 1) * 3) bad++;
    });
  }
  eq('环状 ORF 长度不超一圈且无重复', bad + dup, 0, '最长 ' + maxLen);
}

/* ---------- 4. PCR 产物 ---------- */
console.log('\n=== 4. PCR 产物组合（对拍朴素实现）===');
{
  let cases = 0, mismatch = 0, firstBad = null;
  const fwds = ['ACGTACGTACGT', 'AAAACCCCGGGG', 'GGCATTAGCTAA'];
  for (let seed = 1; seed <= 20; seed++) {
    const t = randSeq(120 + seed * 11, seed * 5171);
    for (const fwd of fwds) {
      const rev = C.revComp(fwd);
      for (const circular of [false, true]) {
        cases++;
        const mine = C.pcrProducts(t, fwd, rev, { circular }).map(p => p.size).sort((a, b) => a - b);
        const want = naiveProducts(t, fwd, rev, circular);
        if (JSON.stringify(mine) !== JSON.stringify(want)) {
          mismatch++;
          if (!firstBad) firstBad = { seed, fwd, circular, mine: mine.slice(0, 6), want: want.slice(0, 6) };
        }
      }
    }
  }
  console.log('  对比 ' + cases + ' 组，不一致 ' + mismatch + ' 组');
  if (firstBad) console.log('  首个不一致：' + JSON.stringify(firstBad));
  eq('PCR 产物对拍', mismatch, 0);
}

/* ---------- 5. 反向互补 / 翻译 逐条校验 ---------- */
console.log('\n=== 5. 反向互补与密码表 ===');
{
  const COMP = { A: 'T', T: 'A', G: 'C', C: 'G' };
  let bad = 0;
  for (let seed = 1; seed <= 30; seed++) {
    const s = randSeq(30 + seed, seed * 101);
    let naive = '';
    for (let i = s.length - 1; i >= 0; i--) naive += COMP[s[i]];
    if (C.revComp(s) !== naive) bad++;
  }
  eq('反向互补 30 组', bad, 0);

  // 密码表：64 个密码子逐一核对
  const STD = {
    TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L', CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
    ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M', GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
    TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S', CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
    ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T', GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
    TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*', CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
    AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K', GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
    TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W', CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
    AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R', GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G'
  };
  let codonBad = [];
  Object.keys(STD).forEach(c => {
    if (C.translate(c) !== STD[c]) codonBad.push(c + '→' + C.translate(c) + ' 应为 ' + STD[c]);
  });
  eq('64 个密码子全部正确', codonBad.length, 0, codonBad.slice(0, 3).join('; '));

  // 大小写与 U/T 混用
  eq('小写输入等价', C.translate('atggccattgtaatgggccgctaa'), C.translate('ATGGCCATTGTAATGGGCCGCTAA'));
  eq('RNA 输入等价', C.translate('AUGGCCAUUGUAAUGGGCCGCUAA'), C.translate('ATGGCCATTGTAATGGGCCGCTAA'));
}

/* ---------- 6. 引物设计的约束必须被真正满足 ---------- */
console.log('\n=== 6. 引物设计输出的约束满足性 ===');
{
  // 用真实质粒序列：随机序列的 Tm/GC 分布很难满足引物约束，得不到有效验证
  const tpl = EX[0].seq;
  const cfg = { minLen: 18, maxLen: 26, minTm: 55, maxTm: 66, minGc: 40, maxGc: 62, maxProduct: 4000, topN: 8 };
  const pairs = C.designPrimers(tpl, cfg);
  let viol = [];
  pairs.forEach((p, i) => {
    [p.fwd, p.rev].forEach((q, k) => {
      if (q.seq.length < cfg.minLen || q.seq.length > cfg.maxLen) viol.push('#' + i + ' 长度 ' + q.seq.length);
      if (q.gc < cfg.minGc - 0.01 || q.gc > cfg.maxGc + 0.01) viol.push('#' + i + ' GC ' + q.gc.toFixed(1));
      if (q.tm < cfg.minTm - 0.01 || q.tm > cfg.maxTm + 0.01) viol.push('#' + i + ' Tm ' + q.tm.toFixed(1));
    });
    if (p.size > cfg.maxProduct) viol.push('#' + i + ' 产物 ' + p.size);
    if (p.tmDiff > 3.0001) viol.push('#' + i + ' Tm 差 ' + p.tmDiff.toFixed(2));
    // 引物必须真的来自模板
    if (tpl.indexOf(p.fwd.seq) === -1) viol.push('#' + i + ' 正向引物不在模板中');
    if (tpl.indexOf(C.revComp(p.rev.seq)) === -1) viol.push('#' + i + ' 反向位点不在模板中');
  });
  console.log('  设计 ' + pairs.length + ' 对，违反约束 ' + viol.length + ' 处');
  if (viol.length) console.log('  ' + viol.slice(0, 5).join('; '));
  eq('引物设计满足全部约束', viol.length, 0);
}

/* ---------- 7. 比对的对拍 ---------- */
console.log('\n=== 7. 双序列比对 ===');
{
  // 相同序列
  const s = randSeq(80, 777);
  eq('相同序列一致度 100%', C.alignPair(s, s, {}).identity, 100);
  eq('相同序列 gaps=0', C.alignPair(s, s, {}).gaps, 0);
  // 已知编辑：删掉中间 10 nt，应产生 10 个 gap
  const del = s.slice(0, 30) + s.slice(40);
  const r = C.alignPair(s, del, {});
  eq('删除 10 nt → 10 个 gap', r.gaps, 10);
  eq('去 gap 后还原', r.alignedB.replace(/-/g, ''), del);
  eq('去 gap 后还原(参考)', r.alignedA.replace(/-/g, ''), s);
  // 三行等长
  eq('三行等长', r.alignedA.length === r.alignedB.length && r.alignedA.length === r.midline.length, true);
  // 得分应为最优：换一种参数再验
  const r2 = C.alignPair(s, del, { match: 1, mismatch: -1, gap: -1 });
  eq('自定义打分下 gap 数不劣化', r2.gaps <= 10, true, 'gaps=' + r2.gaps);
}

console.log('\n' + '='.repeat(60));
if (fails.length) {
  console.log('对拍发现 ' + fails.length + ' 处不一致：\n');
  fails.forEach(f => console.log('  ✘ ' + f.name + '\n     实得: ' + f.got + '\n     应为: ' + f.want + (f.extra ? '\n     备注: ' + f.extra : '')));
  process.exit(1);
} else {
  console.log('全部交叉验证通过 ✔（酶切 / ORF / 环状 / PCR / 翻译 / 引物约束 / 比对）');
  process.exit(0);
}
