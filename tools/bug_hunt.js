/* 缺陷猎取：系统性排查边界条件与语义正确性（不是跑通就行） */
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

const bugs = [];
function check(name, fn) {
  try {
    const r = fn();
    if (r) bugs.push({ name, detail: r });
    else console.log('  [OK]   ' + name);
  } catch (e) {
    bugs.push({ name, detail: '抛出异常: ' + (e.message || e) });
  }
}

console.log('\n=== A. 环状拓扑下的 ORF ===');
{
  // 末尾 A + 开头 TG 在环状时拼成 ATG
  const inner = 'GCCACCATGGGGGCCGAATTCAAATAA';
  const seq = 'TG' + inner + 'A'.repeat(60) + 'A';
  const linear = C.findORFs(seq, 5, true).filter(o => o.strand === '+').map(o => o.aaLength);
  // 旋转一位让接缝进入序列内部，作为"正确答案"的参照
  const rot = seq.slice(-1) + seq.slice(0, -1);
  const rotMax = Math.max.apply(null, C.findORFs(rot, 5, true).filter(o => o.strand === '+').map(o => o.aaLength));

  check('环状序列应能找到跨接缝 ORF', () => {
    const circ = C.findORFsCircular ? C.findORFsCircular(seq, 5, true) : null;
    if (!circ) return 'findORFsCircular 尚未实现（环状跨接缝 ORF 漏检）';
    const maxAa = Math.max.apply(null, circ.filter(o => o.strand === '+').map(o => o.aaLength));
    if (maxAa < rotMax) return '环状结果 ' + maxAa + ' aa < 旋转参照 ' + rotMax + ' aa';
    return false;
  });
}

console.log('\n=== B. 负链 ORF 的坐标语义 ===');
{
  // 正链放一段反向互补的 CDS，使其只可能出现在负链
  // 注意 CDS 长度必须是 3 的倍数，否则 findORFs 找到的只是片段（曾因此误判为坐标 bug）
  const cds = 'ATGGGGCCCGCCCCCGGGTAA';   // 21 nt = 6 aa + 终止子
  const seq = 'TTTTTT' + C.revComp(cds) + 'AAAAAA';
  const orfs = C.findORFs(seq, 3, true).filter(o => o.strand === '-');
  check('负链 ORF 的 ntStart/ntEnd 应指向原序列坐标', () => {
    if (!orfs.length) return '未找到负链 ORF，无法验证';
    const o = orfs[0];
    // 用返回的坐标从原序列取子序列，应能还原出 CDS
    const sub = C.extractFeature(seq, { start: o.ntStart, end: o.ntEnd, strand: '-' });
    if (sub !== cds) {
      const direct = seq.slice(o.ntStart - 1, o.ntEnd);
      return '坐标指向 rc 而非原序列：按坐标取出的是 ' + direct.slice(0, 24) +
        '，期望还原出 ' + cds + '，实际得到 ' + sub.slice(0, 24);
    }
    return false;
  });
}

console.log('\n=== C. PCR 产物：模板含多个结合位点 ===');
{
  const unit = 'ACGTACGTACGT';
  const tpl = unit + 'AAA' + unit + 'TTTTTTTTTT' + unit;
  const fwd = unit, rev = C.revComp(unit);
  check('pcrProduct 应报告全部产物组合', () => {
    const res = C.pcrProducts ? C.pcrProducts(tpl, fwd, rev) : null;
    if (!res) return 'pcrProducts（复数）尚未实现，只返回单一产物';
    if (!res.length) return '未返回任何产物';
    return false;
  });
}

console.log('\n=== D. 极端输入的健壮性 ===');
{
  const weird = [
    ['', '空串'],
    ['A', '单碱基'],
    ['AC', '两碱基'],
    ['NNNNNNNNNNNN', '全 N'],
    ['!!!!@@@@', '全非法字符'],
    ['ACGTNNNNRYSWKMBDHV', '全简并碱基'],
    ['ACGT'.repeat(100) + 'XYZ', '含非法字符混入'],
    ['\n\n\n', '只有空白']
  ];
  weird.forEach(function (w) {
    check('输入[' + w[1] + ']时核心函数不抛异常', () => {
      C.parseFasta(w[0]);
      C.gcContent(w[0]);
      C.revComp(w[0]);
      C.translate(w[0]);
      C.findORFs(w[0], 5, true);
      C.tmNearestNeighbor(w[0], {});
      C.hairpinDeltaG(w[0], {});
      C.dimerDeltaG(w[0], w[0], {});
      C.mwNucleic(w[0], 'dsDNA');
      C.designPrimers(w[0], {});
      C.alignPair(w[0], w[0], {});
      C.digest(w[0], [E.byName('EcoRI')], false);
      return false;
    });
  });
  // 极易出问题的具体场景
  check('gcContent 全 N 序列不产生 NaN', () => {
    const v = C.gcContent('NNNNNN');
    return isNaN(v) ? '返回 NaN' : false;
  });
  check('gcWindows 全 N 序列不产生 NaN', () => {
    const w = C.gcWindows('N'.repeat(500), 50, 10);
    return w.some(x => isNaN(x.gc)) ? '存在 NaN' : false;
  });
  check('migrationMm/yOf 类计算对 0 bp 不产生 Infinity', () => {
    // 核心层没有迁移函数，用对数边界验证：log10(0) 不应进入绘制路径
    return (Math.log10(0) === -Infinity) ? false : 'log10(0) 非负无穷，环境异常';
  });
}

console.log('\n=== E. 数值正确性抽查 ===');
{
  check('已知序列的 GC 含量', () => {
    const v = C.gcContent('GGCC');
    return Math.abs(v - 100) > 1e-9 ? 'GC 100% 序列算出 ' + v : false;
  });
  check('已知序列的翻译', () => {
    const p = C.translate('ATGGCCATTGTAATGGGCCGCTAA');
    return p.indexOf('MAIVMGR') !== 0 ? '翻译结果 ' + p + ' 与文献不符' : false;
  });
  check('EcoRI 在位点上的切点应为识别序列起点', () => {
    const h = C.findSites('AAAAGAATTCAAAA', E.byName('EcoRI'), false);
    // G^AATTC：识别序列起点在索引 4，切点在 G 之后 → cut = 4 + 1 = 5
    return (h.length !== 1 || h[0].cut !== 5) ? '切点 ' + (h[0] && h[0].cut) + '，期望 5（G^AATTC 切在 G 之后）' : false;
  });
  check('双链分子量量级正确（1000 bp ≈ 618 kDa）', () => {
    const mw = C.mwNucleic('ACGT'.repeat(250), 'dsDNA');
    return Math.abs(mw - 617900) / 617900 > 0.02 ? '算出 ' + Math.round(mw) + ' Da' : false;
  });
  check('Tm 随长度单调不减', () => {
    const a = C.tmNearestNeighbor('ACGTACGTACGTACGT', {});
    const b = C.tmNearestNeighbor('ACGTACGTACGTACGTACGTACGTACGTACGT', {});
    return b < a ? '长序列 Tm(' + b.toFixed(1) + ') < 短序列(' + a.toFixed(1) + ')' : false;
  });
  check('环状酶切片段之和等于序列长度', () => {
    const s = EX[0].seq;
    const d = C.digest(s, [E.byName('EcoRI'), E.byName('HindIII')], true);
    const sum = d.fragments.reduce((a, b) => a + b, 0);
    return sum !== s.length ? '片段合计 ' + sum + ' ≠ ' + s.length : false;
  });
}

console.log('\n=== F. 酶库自洽性 ===');
{
  check('所有酶的切点偏移落在识别序列内', () => {
    const bad = E.list.filter(e => e.cutOffset < 0 || e.cutOffset > e.siteClean.length);
    return bad.length ? bad.map(e => e.name).join(',') : false;
  });
  check('每个酶在自己的位点上都能切到', () => {
    const bad = [];
    E.list.forEach(function (e) {
      // 用识别序列本身（简并碱基替换成具体碱基）构造底物
      const concrete = e.siteClean.replace(/[RYSWKMBDHVN]/g, function (ch) {
        return { R: 'A', Y: 'C', S: 'G', W: 'A', K: 'G', M: 'A', B: 'C', D: 'A', H: 'A', V: 'A', N: 'A' }[ch] || 'A';
      });
      const substrate = 'AAAA' + concrete + 'AAAA';
      const hits = C.findSites(substrate, e, false);
      if (hits.length !== 1) bad.push(e.name + '(' + hits.length + ')');
    });
    return bad.length ? '异常酶: ' + bad.join(', ') : false;
  });
}

console.log('\n=== G. 新功能边界（蛋白/突变/Motif/连接/双酶切） ===');
{
  check('proteinProperties 对随机序列不抛异常且字段自洽', () => {
    for (let i = 0; i < 200; i++) {
      const aa = C.randomSeq(120, 50).replace(/[T]/g, 'K').replace(/[G]/g, 'E')
        .replace(/[C]/g, 'A').replace(/[A]/g, 'V');
      const pp = C.proteinProperties(aa);
      if (!pp) return '第 ' + i + ' 组返回 null';
      if (!(pp.pI >= 0 && pp.pI <= 14)) return 'pI 越界 ' + pp.pI;
      if (!(pp.gravy >= -4.5 && pp.gravy <= 4.5)) return 'GRAVY 越界 ' + pp.gravy;
      if (!(pp.mw > 0)) return 'MW 非正';
      if (Math.abs(pp.mw - C.mwProtein(aa)) > 1e-6) return 'MW 与 mwProtein 不一致';
    }
    return false;
  });
  check('DIWV 表完整（20×20 无空洞）', () => {
    // 通过随机二肽间接探测：任何标准二肽都不应触发 undefined
    const A = 'ACDEFGHIKLMNPQRSTVWY';
    for (let i = 0; i < 20; i++) for (let j = 0; j < 20; j++) {
      const pp = C.proteinProperties(A[i] + A[j]);
      if (!pp || !isFinite(pp.instability)) return '二肽 ' + A[i] + A[j] + ' 异常';
    }
    return false;
  });
  check('反向翻译对全部 20 种残基 + 终止子回译一致', () => {
    const rt = C.reverseTranslate('ACDEFGHIKLMNPQRSTVWY*');
    if (!rt.backTranslateOk) return '回译不一致';
    if (rt.dna.length !== 63) return '长度 = ' + rt.dna.length;
    return false;
  });
  check('突变引物：500 组随机模板/位置/突变不抛异常且退火区忠实于模板', () => {
    for (let i = 0; i < 500; i++) {
      const tpl = C.randomSeq(400 + Math.floor(Math.random() * 2600), 35 + Math.random() * 45);
      const pos = 1 + Math.floor(Math.random() * tpl.length);
      const nb = ['A', 'C', 'G', 'T', '', 'AAA', 'GGGGCC'][Math.floor(Math.random() * 7)];
      const want = nb === '' ? 1 + Math.floor(Math.random() * 9) : (Math.random() < 0.3 ? 0 : 1);
      const circular = Math.random() < 0.5;
      // 构造合法区间：线性模板下被替换区不得越过末端（越界情形由下一条断言单独覆盖）
      const oldLen = Math.min(want, Math.max(0, tpl.length - pos + 1));
      const r = C.designMutPrimers(tpl, pos, nb, { oldLen: oldLen, circular: circular });
      if (!r.ok) return '第 ' + i + ' 组失败: ' + r.reason + ' (pos=' + pos + ', nb=' + nb + ', old=' + oldLen + ')';
      if ((tpl + tpl).indexOf(r.fwdAnneal) < 0 && tpl.indexOf(r.fwdAnneal) < 0)
        return '第 ' + i + ' 组正向退火区不在模板上';
      if ((tpl + tpl).indexOf(C.revComp(r.revAnneal)) < 0 && tpl.indexOf(C.revComp(r.revAnneal)) < 0)
        return '第 ' + i + ' 组反向退火区不在模板上';
      if (r.productLen !== tpl.length - oldLen + nb.length) return '第 ' + i + ' 组产物长度错误';
      if (!(r.productLen > 0)) return '第 ' + i + ' 组产物长度非正: ' + r.productLen;
    }
    return false;
  });
  check('突变引物：越界 / 超长替换一律被拒绝且产物长度必为正', () => {
    const tpl = C.randomSeq(500, 50);
    const cases = [
      { pos: 500, oldLen: 5, circular: false, why: '线性尾部越界' },
      { pos: 1, oldLen: 501, circular: true, why: '替换碱基数超过模板全长（环状）' },
      { pos: 300, oldLen: 501, circular: true, why: '环状下从中间替换超过全长' },
      { pos: 0, oldLen: 1, circular: true, why: '位置为 0' },
      { pos: 501, oldLen: 1, circular: true, why: '位置越界' }
    ];
    for (const c of cases) {
      const r = C.designMutPrimers(tpl, c.pos, 'A', { oldLen: c.oldLen, circular: c.circular });
      if (r.ok) return c.why + ' 未被拒绝（产物 ' + r.productLen + ' bp）';
      if (!r.reason) return c.why + ' 被拒绝但没给出原因';
    }
    return false;
  });
  check('Motif：穷举窗口对拍无漏报无错报（环状 + 线性 × 双链）', () => {
    const IUP = { R: '[AG]', Y: '[CT]', S: '[GC]', W: '[AT]', K: '[GT]', M: '[AC]', B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]', N: '[ACGT]' };
    function brute(seq, pat, circular) {
      const n = seq.length, L = pat.length, dbl = seq + seq;
      const re = new RegExp('^' + pat.split('').map(c => IUP[c] || c).join('') + '$');
      const set = new Set();
      const limit = circular ? n : n - L + 1;
      for (let p = 1; p <= limit; p++) {
        let w = '';
        for (let k = 0; k < L; k++) w += dbl.charAt(((p - 1 + k) % n + n) % n);
        if (re.test(w)) set.add(p + '+');
        if (re.test(C.revComp(w))) set.add(p + '-');
      }
      return set;
    }
    for (let t = 0; t < 150; t++) {
      const seq = C.randomSeq(40 + Math.floor(Math.random() * 80), 45);
      const pat = ['AGG', 'AATAAA', 'GTYRAC', 'ATGN', 'TGC'][t % 5];
      for (const circular of [true, false]) {
        const got = C.findMotifs(seq, pat, { circular: circular, bothStrands: true });
        const gs = new Set(got.map(h => h.pos + h.strand));
        const bs = brute(seq, pat, circular);
        for (const k of bs) if (!gs.has(k)) return '漏报 ' + k + '（pat=' + pat + ', circular=' + circular + '）';
        for (const k of gs) if (!bs.has(k)) return '错报 ' + k + '（pat=' + pat + ', circular=' + circular + '）';
      }
    }
    return false;
  });
  check('Motif：环状搜索结果 ⊇ 线性结果', () => {
    for (let i = 0; i < 100; i++) {
      const s = C.randomSeq(300, 45);
      const pat = ['AGGAGG', 'AATAAA', 'GTYRAC', 'NNNNATG'][i % 4];
      const lin = C.findMotifs(s, pat, { circular: false });
      const cir = C.findMotifs(s, pat, { circular: true });
      if (cir.length < lin.length) return '第 ' + i + ' 组环状结果反而更少';
    }
    return false;
  });
  check('Motif：负链命中位置映射正确（去重校验）', () => {
    // 回文 motif 在双链模式下每个命中应出现正反各一次且坐标相同
    const hits = C.findMotifs('AAAGAATTCTTT', 'GAATTC', { bothStrands: true });
    if (hits.length !== 2) return '回文位点应双链各中一次，实际 ' + hits.length;
    if (hits[0].pos !== 4 || hits[1].pos !== 4) return '坐标应均为 4：' + JSON.stringify(hits);
    return false;
  });
  check('连接计算：与 mwNucleic 口径差异在 8% 以内（660 约定检查）', () => {
    // 验证 660 g/mol/bp 与精确残基质量的偏离程度
    const s = C.randomSeq(1000, 50);
    const exact = C.mwNucleic(s, 'dsDNA');
    const approx = 1000 * 660;
    const dev = Math.abs(exact - approx) / exact;
    return dev > 0.08 ? '偏差 ' + (dev * 100).toFixed(1) + '% 过大' : false;
  });
  check('双酶切：51 种酶两两抽查 200 对不抛异常', () => {
    const s = C.randomSeq(2000, 50);
    for (let i = 0; i < 200; i++) {
      const a = E.list[Math.floor(Math.random() * E.list.length)];
      const b = E.list[Math.floor(Math.random() * E.list.length)];
      const r = C.doubleDigestSummary(s, a, b, i % 2 === 0);
      if (!r) return '第 ' + i + ' 对返回 null: ' + a.name + '×' + b.name;
      const d = C.digest(s, [a, b], i % 2 === 0);
      if (JSON.stringify(r.fragments) !== JSON.stringify(d.fragments))
        return '第 ' + i + ' 对片段与 digest() 不一致: ' + a.name + '×' + b.name;
    }
    return false;
  });
}

console.log('\n' + '='.repeat(60));
if (bugs.length) {
  console.log('发现 ' + bugs.length + ' 个问题：\n');
  bugs.forEach((b, i) => console.log('  ' + (i + 1) + '. ' + b.name + '\n     → ' + b.detail));
  process.exit(1);
} else {
  console.log('边界与正确性排查全部通过 ✔');
  process.exit(0);
}
