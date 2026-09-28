/* MolBench 冒烟测试：在 node 中加载浏览器脚本并验证核心算法 */
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const root = path.join(__dirname, '..');

function loadBrowser(f) {
  const ctx = { window: {}, console };
  vm.createContext(ctx);
  vm.runInContext(fs.readFileSync(path.join(root, 'assets', f), 'utf8'), ctx, { filename: f });
  return ctx.window;
}

const C = loadBrowser('core.js').MolCore;
const E = loadBrowser('enzymes.js').MolEnzymes;
const EX = loadBrowser('examples.js').MolExamples;

let fails = 0;
function ok(cond, label, extra) {
  console.log((cond ? '  [OK]   ' : '  [FAIL] ') + label + (extra !== undefined ? '  → ' + extra : ''));
  if (!cond) fails++;
}
function near(a, b, tol, label) {
  ok(Math.abs(a - b) <= tol, label, a.toFixed(3) + ' ≈ ' + b);
}

console.log('\n--- 1. 序列解析 ---');
const fa = '>test gene\n' + 'ACGT\n'.repeat(5);
const p = C.parseFasta(fa);
ok(p.name === 'test gene', 'FASTA 名称解析', p.name);
ok(p.seq === 'ACGT'.repeat(5), 'FASTA 序列拼接', p.seq.length + ' nt');
ok(C.detectType('ACGGTA') === 'dna', 'DNA 识别');
ok(C.detectType('ACGGUA') === 'rna', 'RNA 识别');
ok(C.detectType('MKWVTFISLLLLFSSAYS') === 'protein', '蛋白识别');

// GenBank / EMBL：只取 ORIGIN / SQ 段，注释里的英文不得混入序列
const GB = [
  'LOCUS       TEST        40 bp    DNA     linear   UNK 01-JAN-2026',
  'DEFINITION  test sequence for unit test.',
  'ACCESSION   TEST0001',
  'FEATURES             Location/Qualifiers',
  '     CDS             1..40',
  '                     /gene="test"',
  '                     /product="test protein"',
  'ORIGIN',
  '        1 atggccattg taatgggccg ctgaaagggt gcccgatagc',
  '//'
].join('\n');
const gb = C.parseFasta(GB);
ok(gb.seq.length === 40, 'GenBank 只解析 ORIGIN 段（注释不计入）', gb.seq.length + ' nt');
ok(gb.seq === 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGC', 'GenBank 序列内容正确', gb.seq);
ok(gb.name === 'TEST', 'GenBank 名称取自 LOCUS', gb.name);
ok(gb.type === 'dna', 'GenBank 类型识别');

const EMBL = [
  'ID   TEST; SV 1; linear; DNA; STD; UNC; 40 BP.',
  'DE   test embl sequence',
  'SQ   Sequence 40 BP;',
  '     atggccattg taatgggccg ctgaaagggt gcccgatagc',
  '//'
].join('\n');
const em = C.parseFasta(EMBL);
ok(em.seq.length === 40, 'EMBL 只解析 SQ 段', em.seq.length + ' nt');
ok(em.seq === 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGC', 'EMBL 序列内容正确');

// FEATURES 段解析
const GBF = [
  'LOCUS       TEST       100 bp    DNA     linear   UNK 01-JAN-2026',
  'FEATURES             Location/Qualifiers',
  '     source          1..100',
  '                     /organism="Escherichia coli"',
  '     gene            10..90',
  '                     /gene="testGene"',
  '     CDS             complement(20..88)',
  '                     /gene="testGene"',
  '                     /product="test protein"',
  '                     /translation="MKWVTFISLLLLFSSAYSRGVF',
  '                     RRDTHKSEIAHRFKDLGEEH"',
  '     mRNA            join(1..30,50..100)',
  '                     /gene="testGene"',
  'ORIGIN',
  '        1 acgtacgtac',
  '//'
].join('\n');
const gbf = C.parseFasta(GBF).features;
ok(gbf.length === 4, 'FEATURES 解析出 4 条注释', gbf.length + ' 条');
const cdsF = gbf.filter(f => f.type === 'CDS')[0];
ok(!!cdsF, '找到 CDS 注释');
ok(cdsF && cdsF.strand === '-', 'complement 识别为负链', cdsF && cdsF.strand);
ok(cdsF && cdsF.start === 20 && cdsF.end === 88, 'complement 区间解析', cdsF && (cdsF.start + '-' + cdsF.end));
ok(cdsF && cdsF.qualifiers.translation.length === 42, '跨多行 /translation 正确拼接', cdsF && cdsF.qualifiers.translation.length + ' aa');
ok(cdsF && cdsF.qualifiers.product === 'test protein', 'product qualifier 解析');
const mrnaF = gbf.filter(f => f.type === 'mRNA')[0];
ok(mrnaF && mrnaF.join === true, 'join 位置被标记', mrnaF && String(mrnaF.join));
ok(mrnaF && mrnaF.start === 1 && mrnaF.end === 100, 'join 区间取首尾', mrnaF && (mrnaF.start + '-' + mrnaF.end));
// extractFeature：负链取反向互补
const gseq = C.parseFasta(GBF).seq;
const subF = C.extractFeature('ACGTACGTACGTACGTACGTACGTACGTACGT', { start: 1, end: 4, strand: '-' });
ok(subF === 'ACGT', 'extractFeature 正链', subF);
ok(C.extractFeature('ACGTACGT', { start: 1, end: 4, strand: '-' }) === 'ACGT', 'extractFeature 负链取自反向互补');
ok(C.extractFeature('ACGTACGT', { start: 100, end: 200, strand: '+' }) === '', '区间越界返回空串');

// 格式识别不互相误判
ok(C.detectFormat('>a\nACGT') === 'fasta', '格式识别 fasta');
ok(C.detectFormat(GB) === 'genbank', '格式识别 genbank');
ok(C.detectFormat(EMBL) === 'embl', '格式识别 embl');
ok(C.detectFormat('ACGTACGT') === 'raw', '格式识别 raw');
// 没有 ORIGIN 段的 GenBank 不应崩溃（回退通用清洗）
ok(typeof C.parseFasta('LOCUS x\nDEFINITION nothing here.').seq === 'string', '无 ORIGIN 段时不崩溃');

console.log('\n--- 2. 反向互补 / 翻译 ---');
ok(C.revComp('ATGC') === 'GCAT', 'revComp ATGC → GCAT', C.revComp('ATGC'));
ok(C.revComp('GAATTC') === 'GAATTC', 'EcoRI 位点回文自反');
ok(C.translate('ATGGCCATTGTAATGGGCCGCTAA') === 'MAIVMGR*', '标准翻译 (文献用例)', C.translate('ATGGCCATTGTAATGGGCCGCTAA'));
ok(C.transcribe('ATGC') === 'AUGC', '转录 T→U');

console.log('\n--- 3. GC / Tm ---');
near(C.gcContent('GGCC'), 100, 0.001, 'GC 100%');
near(C.gcContent('ATAT'), 0, 0.001, 'GC 0%');
near(C.tmWallace('ATGC'), 12, 0.001, 'Wallace ATGC = 2*2+4*2 = 12');
// 文献对照：Primer3 对 GTAAAACGACGGCCAGT (M13 正向引物, 17nt, 50% GC) Tm ≈ 55.5~58
// M13 通用引物，IDT OligoAnalyzer 参考值约 52–56 °C
const m13 = 'GTAAAACGACGGCCAGT';
const tmNN = C.tmNearestNeighbor(m13, { conc: 500, na: 50, mg: 1.5, dntp: 0.2 });
ok(tmNN > 48 && tmNN < 62, 'M13 正向引物 Tm 落在合理区间', tmNN.toFixed(2) + ' °C');
const tm47 = C.tmNearestNeighbor('CAGGAAACAGCTATGAC', {});
ok(tm47 > 48 && tm47 < 62, 'M13-47 反向引物 Tm 落在合理区间', tm47.toFixed(2) + ' °C');
// 长序列 Tm 应高于短序列；GC 高者 Tm 更高
ok(C.tmNearestNeighbor('GCGCGCGCGCGCGCGCGCGC', {}) > C.tmNearestNeighbor('ATATATATATATATATATAT', {}),
   'GC 富集引物 Tm 高于 AT 富集引物');

console.log('\n--- 4. 分子量 ---');
// 1000 bp 随机序列 dsDNA MW 应 ≈ 618,000
const r1k = 'ACGT'.repeat(250);
const mwds = C.mwNucleic(r1k, 'dsDNA');
ok(Math.abs(mwds - 617900) / 617900 < 0.02, 'dsDNA 1000bp MW ≈ 618 kDa', (mwds / 1000).toFixed(1) + ' kDa');
// 单链 24nt 引物
const mwss = C.mwNucleic(m13, 'ssDNA');
ok(mwss > 5000 && mwss < 8000, 'ssDNA 17nt MW 合理', mwss.toFixed(1) + ' Da');
// 蛋白：MAIVMGR 七肽，手算 = 131.19+71.08+113.16+99.13+131.19+57.05+156.19+18.02 = 777.0
const mwp = C.mwProtein('MAIVMGR');
near(mwp, 777.0, 0.5, '多肽 MAIVMGR MW 与手算一致');

console.log('\n--- 5. 浓度换算 ---');
const mw1kb = 617900;
near(C.concNgUlToNm(C.concNmToNgUl(100, mw1kb), mw1kb), 100, 0.01, 'nM ↔ ng/µL 往返一致');
// 1000bp dsDNA 1 nM = 0.6179 ng/µL
near(C.concNmToNgUl(1, mw1kb), 0.6179, 0.001, '1 nM 1kb dsDNA = 0.6179 ng/µL');

console.log('\n--- 6. 酶切 ---');
const seq1 = 'AAAAGAATTCAAAAAGCTTTT';
const eco = E.byName('EcoRI'), hin = E.byName('HindIII');
const d1 = C.digest(seq1, [eco], false);
ok(d1.perEnzyme[0].count === 1, '线性序列中 EcoRI 切点数', d1.perEnzyme[0].count);
ok(d1.fragments.length === 2, '单切点 → 2 片段', d1.fragments.join('+'));
ok(d1.fragments.reduce((a, b) => a + b, 0) === seq1.length, '片段总和 = 原长');
const d2 = C.digest(seq1, [eco, hin], false);
ok(d2.fragments.length === 3, '双酶切 → 3 片段', d2.fragments.join('+'));
// 环状单切点 → 1 个全长片段
const circ = 'GAATTC' + 'A'.repeat(94);
const d3 = C.digest(circ, [eco], true);
ok(d3.fragments.length === 1 && d3.fragments[0] === 100, '环状单切 → 线性化 100 bp', d3.fragments.join('+'));
// 环状双切点 → 2 片段
const circ2 = 'GAATTC' + 'A'.repeat(40) + 'GAATTC' + 'A'.repeat(48);
const d4 = C.digest(circ2, [eco], true);
ok(d4.fragments.length === 2 && d4.fragments.reduce((a, b) => a + b, 0) === 100, '环状双切 → 2 片段合计 100', d4.fragments.join('+'));
// 简并酶：AvaI = C^YCGRG，CCCGGG / CTCGAG 均为其位点
const avaI = E.byName('AvaI');
const avaHit = C.findSites('TTTTCCCGGGTTTT', avaI);
ok(avaHit.length === 1 && avaHit[0].matched === 'CCCGGG', '简并酶 AvaI (C^YCGRG) 匹配', JSON.stringify(avaHit.map(h => h.matched)));
ok(C.findSites('AAACTCGAGAAA', avaI).length === 1, '简并酶 AvaI 识别 CTCGAG');
// 简并酶不应误配非位点
ok(C.findSites('AAAGGGCCCAAA', avaI).length === 0, '简并酶不误配 GGGCCC');

console.log('\n--- 7. 示例数据集 ---');
const pl = EX[0].seq;
const orfs = C.findORFs(pl, 60, true);
ok(orfs.length > 0, '质粒中检出 ORF', orfs.length + ' 个');
// makeCDS(148) → 149 个密码子 = 148 aa + 终止子
ok(orfs.some(o => o.aaLength === 148), '检出的 ORF 含设计的 148 aa CDS', orfs.map(o => o.aaLength).join('/'));
ok(orfs.every(o => o.protein.indexOf('*') === -1), '所有 ORF 内部无终止子');
ok(orfs.every(o => o.protein.charAt(0) === 'M'), '所有 ORF 以 M 起始');
ok(orfs.every(o => o.nt.length === (o.aaLength + 1) * 3), 'ORF 核酸长度 = (aa+1)×3（含终止子）');
const dbl = C.digest(pl, [eco, hin], true);
ok(dbl.fragments.reduce((a, b) => a + b, 0) === pl.length, '质粒双酶切片段总和 = 全长', dbl.fragments.join(' + '));
const gw = C.gcWindows(pl, 50, 10);
ok(gw.length > 0 && gw.every(x => x.gc >= 0 && x.gc <= 100), 'GC 窗口计算合法', gw.length + ' 个窗口');

console.log('\n--- 8. PCR 产物 ---');
const tpl = 'ACGTA'.repeat(5) + 'ACGTACGTACGTACGTTACGAAAGCTT' + 'ACGTA'.repeat(20) + 'TTGGCACGTACGTACGT' + 'ACGTA'.repeat(5);
const fwd = 'ACGTA'.repeat(5).slice(0, 18);
const revSite = 'TTGGCACGTACGTACGT';
const rev = C.revComp(revSite);
const prs = C.pcrProducts(tpl, fwd, rev, {});
ok(prs.length > 0, 'PCR 产物定位成功', prs.length ? prs.length + ' 种，最短 ' + prs[0].size + ' bp' : '未找到');
const p0 = prs[0];
const p0seq = C.pcrProductSeq(tpl, p0);
ok(p0seq.slice(0, 18) === fwd, '产物 5′ 端 = 正向引物');
ok(p0seq.slice(-revSite.length) === revSite, '产物 3′ 端 = 反向引物结合位点');
ok(p0seq.length === p0.size, '产物序列长度与报告的大小一致', p0seq.length + ' vs ' + p0.size);

console.log('\n--- 9. 引物风险 ---');
ok(C.gcClamp('ACGTGCGGCC').gc >= 3, 'GC clamp 检测', C.gcClamp('ACGTGCGGCC').gc + '/5');
// 旧的非热力学实现（maxComplementaryRun / threePrimeDimer / hairpin / pcrProduct 单数版）
// 已被 ΔG 版本取代并删除，避免两套实现并存导致漂移；这里统一用 ΔG 版断言
ok(C.dimerDeltaG('GCGCGCGC', 'GCGCGCGC', {}).run >= 6, '自身二聚体检测（回文）', C.dimerDeltaG('GCGCGCGC', 'GCGCGCGC', {}).run + ' bp');
ok(typeof C.dimerDeltaG('ACGTACGTAAA', 'ACGTACGT', {}).dG === 'number', '交叉二聚体 ΔG 可用');

console.log('\n--- 10. GC 窗口增量算法（回归：曾为 O(n×win)） ---');
{
  const seqs = [C.randomSeq(500, 50), C.randomSeq(2000, 35), 'ACGTN'.repeat(300), 'N'.repeat(50) + C.randomSeq(800, 60)];
  let maxDiff = 0, compared = 0;
  seqs.forEach(sq => {
    [[10, 5], [50, 10], [100, 10], [7, 3]].forEach(([w, st]) => {
      const got = C.gcWindows(sq, w, st).map(x => x.gc);
      const want = [];
      for (let i = 0; i + w <= sq.length; i += st) want.push(C.gcContent(sq.substr(i, w)));
      ok(got.length === want.length, `窗口数量一致 win=${w} step=${st}`, got.length);
      for (let i = 0; i < got.length; i++) { maxDiff = Math.max(maxDiff, Math.abs(got[i] - want[i])); compared++; }
    });
  });
  ok(maxDiff < 1e-9, '增量算法与朴素实现数值一致', compared + ' 个窗口，最大偏差 ' + maxDiff.toExponential(2));
}

console.log('\n--- 11. 环状跨接缝酶切位点 ---');
{
  const eco = E.byName('EcoRI'); // GAATTC
  // 序列末尾 ...GAAT + 开头 TC...，环状时接缝处拼成 GAATTC
  const wrapSeq = 'TC' + 'A'.repeat(40) + 'GAAT';
  const w1 = C.findSites(wrapSeq, eco, true);
  ok(w1.length === 1 && w1[0].wrap === true, '环状模式识别跨接缝位点', JSON.stringify(w1.map(h => ({ s: h.start, cut: h.cut, wrap: h.wrap }))));
  ok(C.findSites(wrapSeq, eco, false).length === 0, '线性模式不识别跨接缝位点');
  const dw = C.digest(wrapSeq, [eco], true);
  ok(dw.fragments.length === 1 && dw.fragments[0] === wrapSeq.length, '跨接缝单切点 → 线性化 1 个全长片段', dw.fragments.join('+'));
  // 普通位点不应被误标记为跨接缝
  const normal = 'A'.repeat(10) + 'GAATTC' + 'A'.repeat(30);
  const wn = C.findSites(normal, eco, true);
  ok(wn.length === 1 && wn[0].wrap === false, '普通位点不误标为跨接缝', JSON.stringify(wn.map(h => ({ s: h.start, wrap: h.wrap }))));
  // 4 碱基酶（短位点）跨接缝
  const mbo = E.byName('MboI'); // ^GATC
  const wm = 'ATC' + 'A'.repeat(20) + 'G';
  ok(C.findSites(wm, mbo, true).length === 1, '短位点酶同样支持跨接缝', JSON.stringify(C.findSites(wm, mbo, true).map(h => h.matched)));
}

console.log('\n--- 12. 性能护栏（1 Mb 序列） ---');
{
  const big = C.randomSeq(1000000, 50);
  let t = Date.now(); C.gcWindows(big, 100, 10); const tGc = Date.now() - t;
  t = Date.now(); C.findSites(big, E.byName('EcoRI'), true); const tFind = Date.now() - t;
  ok(tGc < 200, '1 Mb GC 窗口 < 200 ms', tGc + ' ms');
  ok(tFind < 100, '1 Mb 单酶查找 < 100 ms', tFind + ' ms');
}

console.log('\n--- 13. 自动引物设计 ---');
{
  const tpl = EX[0].seq;
  const cfg = { minLen: 18, maxLen: 27, minTm: 55, maxTm: 65, minGc: 40, maxGc: 60, maxProduct: 4000, topN: 5 };
  const t0 = Date.now();
  const pairs = C.designPrimers(tpl, cfg);
  const ms = Date.now() - t0;
  ok(pairs.length > 0, '能设计出候选引物对', pairs.length + ' 对 / ' + ms + ' ms');
  ok(ms < 3000, '设计耗时可接受', ms + ' ms');

  if (pairs.length) {
    const p = pairs[0];
    ok(tpl.indexOf(p.fwd.seq) !== -1, '正向引物序列确实取自模板', p.fwd.seq);
    ok(tpl.indexOf(C.revComp(p.rev.seq)) !== -1, '反向引物结合位点确实在模板上');
    ok(p.fwd.tm >= cfg.minTm && p.fwd.tm <= cfg.maxTm, '正向 Tm 落在设定区间', p.fwd.tm.toFixed(2));
    ok(p.rev.tm >= cfg.minTm && p.rev.tm <= cfg.maxTm, '反向 Tm 落在设定区间', p.rev.tm.toFixed(2));
    ok(p.fwd.gc >= cfg.minGc && p.fwd.gc <= cfg.maxGc, '正向 GC% 落在设定区间', p.fwd.gc.toFixed(1));
    ok(p.tmDiff <= 3, '两引物 Tm 差 ≤ 3 °C', p.tmDiff.toFixed(2));
    // 最关键的闭环：设计结果必须与产物预测一致
    const prs2 = C.pcrProducts(tpl, p.fwd.seq, p.rev.seq, {});
    ok(prs2.length > 0, '设计出的引物对能定位到产物', prs2.length + ' 种');
    ok(prs2.some(function (x) { return x.size === p.size; }), '设计产物大小出现在预测产物中', JSON.stringify(prs2.map(function (x) { return x.size; })) + ' 含 ' + p.size);
    // 评分降序
    let sorted = true;
    for (let i = 1; i < pairs.length; i++) if (pairs[i].score > pairs[i - 1].score + 1e-9) sorted = false;
    ok(sorted, '候选按评分降序返回');
  }

  // 极端参数应优雅返回空而非报错
  ok(C.designPrimers(tpl, { minTm: 95, maxTm: 99 }).length === 0, 'Tm 条件过苛时返回空');
  ok(C.designPrimers('ACGT', cfg).length === 0, '序列过短时返回空');
  ok(C.designPrimers('', cfg).length === 0, '空序列返回空');
  ok(C.designPrimers(tpl, { start: 0, end: 50 }).length === 0, '目标区域过小时返回空');
}

console.log('\n--- 14. 二级结构 ΔG 热力学 ---');
{
  // Tm 拆分后必须与拆分前一致（回归：nnSum 抽取时曾把 tmNearestNeighbor 改坏）
  const m13 = 'GTAAAACGACGGCCAGT';
  near(C.tmNearestNeighbor(m13, {}), 59.05, 0.05, 'Tm 拆分后结果不变');
  ok(C.tmNearestNeighbor('', {}) !== C.tmNearestNeighbor('', {}) || isNaN(C.tmNearestNeighbor('', {})), '空序列 Tm 为 NaN');

  // ΔG 越负越稳定：长序列要比短序列更负
  const dgShort = C.nnDeltaG37('ACGTACGT');
  const dgLong = C.nnDeltaG37('ACGTACGTACGTACGTACGT');
  ok(dgLong < dgShort, 'ΔG 随序列变长而更负', dgShort.toFixed(2) + ' → ' + dgLong.toFixed(2));

  // 发夹：GC 富集序列应检出稳定茎环，AT 富集不应
  const hpGc = C.hairpinDeltaG('GCGCGCGCGCGCGCGCGCGC', {});
  ok(hpGc.dG <= -3, 'GC 富集序列检出稳定发夹', 'ΔG=' + hpGc.dG.toFixed(2) + ', 茎=' + hpGc.stem + 'bp');
  ok(hpGc.stem >= 6, '发夹茎长合理', hpGc.stem + ' bp');
  const hpAt = C.hairpinDeltaG('AAAAATTTTTAAAAA', {});
  ok(hpAt.dG >= -1.5, 'AT 富集序列不误判为强发夹', 'ΔG=' + hpAt.dG.toFixed(2));
  ok(C.hairpinDeltaG('ACGT', {}).dG === 0, '过短序列不发夹（返回 0）');
  ok(C.hairpinDeltaG('', {}).dG === 0, '空序列不发夹（返回 0）');

  // 二聚体：回文序列应检出 3′ 端二聚体
  const dp = C.dimerDeltaG('ACGTACGT', 'ACGTACGT', {});
  ok(dp.dG <= -5, '回文序列检出自身二聚体', 'ΔG=' + dp.dG.toFixed(2) + ', 配对=' + dp.run + 'bp');
  ok(C.dimerTouches3Prime('ACGTACGT', 'ACGTACGT', dp) === true, '二聚体触及 3′ 端被标记');
  // 两条无互补关系的引物不应检出二聚体
  const dpNone = C.dimerDeltaG('AAAAAAAAAAAA', 'AAAAAAAAAAAA', {});
  ok(dpNone.dG === 0, '无互补关系时不报二聚体', 'ΔG=' + dpNone.dG.toFixed(2));
  ok(C.dimerDeltaG('', '', {}).dG === 0, '空序列二聚体 ΔG 为 0');
  ok(C.dimerDeltaG('NNNN', 'NNNN', {}).dG === 0, '全 N 序列不误报');

  // 自洽性：ΔG = ΔH − T·ΔS/1000
  const s1 = 'GTGCTGAGCTAGCTGATCGGATC';
  const sum = C.nnSum(s1, {});
  near(C.nnDeltaG37(s1, {}), sum[0] - 310.15 * sum[1] / 1000, 0.01, 'ΔG 与 ΔH/ΔS 换算自洽');
}

console.log('\n--- 15. 自动引物设计（ΔG 与模板上下文） ---');
{
  const tpl = EX[0].seq;
  const cfg = { minLen: 18, maxLen: 27, minTm: 55, maxTm: 65, minGc: 40, maxGc: 60, maxProduct: 4000, topN: 5 };
  const t0 = Date.now();
  const pairs = C.designPrimers(tpl, cfg);
  const ms = Date.now() - t0;
  ok(pairs.length > 0, '仍能设计出引物对', pairs.length + ' 对 / ' + ms + ' ms');
  ok(ms < 3000, '加入 ΔG 后仍在性能护栏内', ms + ' ms');
  const p = pairs[0];
  ok(typeof p.crossDg === 'number' && p.crossDg <= 0.0001, '交叉二聚体给出 ΔG', p.crossDg.toFixed(2));
  ok(typeof p.fwd.localGc === 'number' && p.fwd.localGc >= 0 && p.fwd.localGc <= 100, '模板局部 GC 已计算', p.fwd.localGc.toFixed(1) + '%');
  ok(typeof p.fwd.hairpinDg === 'number', '候选引物带发夹 ΔG', p.fwd.hairpinDg.toFixed(2));
  ok(typeof p.touch3 === 'boolean', '带 3′ 端二聚体标记');
  const prs2 = C.pcrProducts(tpl, p.fwd.seq, p.rev.seq, {});
  ok(prs2.some(function (x) { return x.size === p.size; }), '设计结果仍与产物预测闭环一致', JSON.stringify(prs2.map(function (x) { return x.size; })) + ' 含 ' + p.size);
  // 高 GC 模板应给出可解释的空结果，而不是崩溃
  const hi = C.randomSeq(2000, 0.85);
  ok(C.designPrimers(hi, cfg).length === 0, '高 GC 模板返回空（GC 过滤生效）');
  ok(C.designPrimers(hi, Object.assign({}, cfg, { maxGc: 90 })).length >= 0, '放宽 GC 上限后不报错');
}

console.log('\n--- 16. 双序列全局比对 ---');
{
  // 完全相同
  const same = C.alignPair('ATGGCCATTGTAATGGGCCGCTAA', 'ATGGCCATTGTAATGGGCCGCTAA');
  near(same.identity, 100, 0.001, '相同序列一致度 100%');
  ok(same.gaps === 0, '相同序列无 gap');
  ok(same.score === 48, '相同序列得分 = 长度 × match', same.score);

  // 缺失：A 比 B 多 5 nt，应产生 5 个 gap 且分值最优
  const del = C.alignPair('ATGGCCATTGTAATGGGCCG', 'ATGGCCATTGGGCCG');
  ok(del.gaps === 5, '缺失 5 nt 完全以 gap 表示', del.gaps + ' 个 gap');
  ok(del.score === 20, '缺失情形的得分为最优解', del.score);
  ok(del.alignedA.replace(/-/g, '') === 'ATGGCCATTGTAATGGGCCG', '比对后 A 去 gap 等于原序列');
  ok(del.alignedB.replace(/-/g, '') === 'ATGGCCATTGGGCCG', '比对后 B 去 gap 等于原序列');
  ok(del.alignedA.length === del.alignedB.length && del.alignedA.length === del.midline.length, '三行等长');

  // 插入与缺失应互为镜像
  const ins = C.alignPair('ATGGCCATTGGGCCG', 'ATGGCCATTGTAATGGGCCG');
  ok(ins.score === del.score && ins.gaps === del.gaps, '插入与缺失结果对称');

  // 完全不同的序列
  const diff = C.alignPair('AAAA', 'TTTT');
  near(diff.identity, 0, 0.001, '无相同碱基时一致度 0%');

  // 边界
  ok(C.alignPair('', 'ACGT') === null, '空序列返回 null');
  ok(C.alignPair('ACGT', '') === null, '另一侧空序列返回 null');
  const single = C.alignPair('A', 'A');
  ok(single.identity === 100 && single.score === 2, '单碱基比对正确');
  // alignPair 是通用比对（不限定核酸），因此保留所有字母，仅过滤数字与符号
  const withJunk = C.alignPair('ACGTNX123-*', 'ACGTNX');
  ok(withJunk.lenA === 6 && withJunk.lenB === 6, '数字与符号被过滤、字母保留', withJunk.lenA + ' / ' + withJunk.lenB);
  ok(withJunk.identity === 100, '过滤后仍正确比对', withJunk.identity + '%');
  // 小写与大写应等价
  ok(C.alignPair('acgt', 'ACGT').identity === 100, '大小写不敏感');

  // 超限拒绝
  const tooBig = C.alignPair(C.randomSeq(3000), C.randomSeq(3000));
  ok(tooBig && tooBig.tooBig === true, '超长序列被拒绝并给出 tooBig 标记', tooBig && (tooBig.cells + ' cells'));
  // 性能护栏
  let t = Date.now(); C.alignPair(C.randomSeq(1000), C.randomSeq(1000)); const ms = Date.now() - t;
  ok(ms < 1500, '1000×1000 比对在护栏内', ms + ' ms');
}

console.log('\n--- 16b. 多序列比对 ---');
{
  // 核心不变量：所有行等长，且去掉 gap 后必须能还原成原始序列
  function invariants(seqs, label) {
    const r = C.alignMultiple(seqs, {});
    if (!r || r.tooBig) return label + '：未返回结果';
    const lens = new Set(r.rows.map(x => x.length));
    if (lens.size !== 1) return label + '：各行长度不一致 ' + [...lens].join('/');
    for (let i = 0; i < seqs.length; i++) {
      if (r.rows[i].replace(/-/g, '') !== seqs[i].toUpperCase()) {
        return label + '：第 ' + (i + 1) + ' 行去 gap 后无法还原（' + r.rows[i] + '）';
      }
    }
    return false;
  }

  ok(!invariants(['ACGTACGT', 'ACGTACGT', 'ACGTACGT'], '相同序列'), 'MSA 相同序列：等长且可还原');
  ok(!invariants(['ACGTACGT', 'ACGAACGT', 'ACGTACGT'], '替换'), 'MSA 含替换：等长且可还原');
  ok(!invariants(['ACGTACGT', 'ACGACGT', 'ACGTACGT'], '缺失'), 'MSA 含缺失：等长且可还原');
  ok(!invariants(['ACGTACGT', 'ACGTAACGT', 'ACGTACGT'], '插入'), 'MSA 含插入：等长且可还原');
  ok(!invariants(['ATGGCCATTGTAATGGGC', 'ATGGCCATTG', 'ATGGCCATTGTAATGGGCCGCTAA'], '长度差异大'),
    'MSA 长度差异大：等长且可还原');
  ok(!invariants(['AAACCCGGG', 'AAACCCGGGTTT', 'AACCCGGG', 'AAACCCGGGTT'], '四条混合'),
    'MSA 四条混合：等长且可还原');

  // 随机压力：插入/缺失/替换/截断各一条，共 4 条
  let viol = 0, trials = 0;
  function rs(n, seed) {
    let s = '', x = seed;
    for (let i = 0; i < n; i++) { x = (x * 1103515245 + 12345) & 0x7fffffff; s += 'ACGT'[x % 4]; }
    return s;
  }
  for (let seed = 1; seed <= 30; seed++) {
    const base = rs(30, seed);
    const seqs = [base, base.slice(0, 10) + base.slice(11), base.slice(0, 15) + 'A' + base.slice(15),
      base.slice(0, 20) + 'T' + base.slice(21), base.slice(5)];
    trials++;
    if (invariants(seqs, 'seed' + seed)) viol++;
  }
  ok(viol === 0, 'MSA 随机压力测试（' + trials + ' 组）不变量零违规');

  // 保守列与一致度
  const same = C.alignMultiple(['ACGTACGT', 'ACGTACGT', 'ACGTACGT'], {});
  ok(same.conserved === 8 && same.columns === 8, '全同序列：保守列 = 列数', same.conserved + '/' + same.columns);
  near(same.consensusPct, 100, 0.01, '全同序列一致度 100%');
  const diff = C.alignMultiple(['ACGTACGT', 'ACGAACGT', 'ACGTACGT'], {});
  ok(diff.conserved === 7, '单点替换：保守列 = 7', diff.conserved);
  ok(diff.count === 3, '序列条数正确');
  ok(typeof diff.center === 'number' && diff.center >= 0 && diff.center < 3, '中心索引落在范围内', diff.center);
  ok(diff.pairwise.length === 3, '三条序列产生 3 组两两比对', diff.pairwise.length);

  // 边界
  ok(C.alignMultiple(['ACGT'], {}) === null, '单条序列返回 null');
  ok(C.alignMultiple([], {}) === null, '空数组返回 null');
  ok(C.alignMultiple(['', ''], {}) === null, '全空串返回 null');
  const withEmpty = C.alignMultiple(['ACGT', '', 'ACGT'], {});
  ok(withEmpty.count === 2, '空串被过滤后按 2 条处理', withEmpty.count);
  const tooMany = C.alignMultiple(new Array(25).fill('ACGTACGT'), {});
  ok(tooMany && tooMany.tooBig === true, '超过条数上限时拒绝', tooMany && tooMany.reason);
  // 超长序列应被矩阵上限挡住
  const big = C.alignMultiple([rs(3000, 1), rs(3000, 2), rs(3000, 3)], {});
  ok(big && big.tooBig === true, '超长序列被矩阵上限拒绝');
  // 小写输入应被规范化
  ok(!invariants(['acgtacgt', 'ACGTACGT', 'acgtacgt'], '大小写混合'), 'MSA 大小写混合可正常工作');
}

console.log('\n--- 17. GenBank location 解析 ---');
{
  const cases = [
    ['123..456', 1, '+'],
    ['87', 1, '+'],
    ['complement(20..88)', 1, '-'],
    ['join(1..30,50..100)', 2, '+'],
    ['complement(join(10..20,40..50))', 2, '-'],
    ['join(1..10,20..30,40..50)', 3, '+'],
    ['<1..>100', 1, '+'],
    ['1..>50', 1, '+']
  ];
  cases.forEach(([loc, nseg, strand]) => {
    const r = C.parseLocation(loc);
    ok(r.segments.length === nseg && r.strand === strand,
      'location 解析 ' + loc, r.segments.length + ' 段 / 链 ' + r.strand);
  });
  const j = C.parseLocation('join(1..10,40..49)');
  ok(j.join === true, 'join 被标记');
  ok(j.segments[0].start === 1 && j.segments[1].end === 49, 'join 分段边界正确');
  // 尖括号不应把区间拆成两个单碱基段
  const angle = C.parseLocation('<1..>100');
  ok(angle.segments.length === 1 && angle.segments[0].start === 1 && angle.segments[0].end === 100,
    '尖括号端点不拆段', JSON.stringify(angle.segments));

  // join 提取：内含子不计入
  const seq = 'AAAAAAAAAA' + 'C'.repeat(29) + 'GGGGGGGGGG';
  const f = C.parseLocation('join(1..10,40..49)');
  const got = C.extractFeature(seq, { segments: f.segments, strand: f.strand, start: 1, end: 49 });
  ok(got === 'AAAAAAAAAAGGGGGGGGGG', 'join 提取跳过内含子', got.length + ' nt');
  const fc = C.parseLocation('complement(join(1..10,40..49))');
  const gotc = C.extractFeature(seq, { segments: fc.segments, strand: fc.strand, start: 1, end: 49 });
  ok(gotc === C.revComp(got), '负链 join = 正链结果的反向互补');

  // feature 的 length 应为各段之和，而不是首尾跨度
  const GBJ = [
    'LOCUS       T        100 bp',
    'FEATURES             Location/Qualifiers',
    '     CDS             join(1..10,40..49)',
    '                     /gene="x"',
    'ORIGIN',
    '        1 acgtacgtac',
    '//'
  ].join('\n');
  const cdsJ = C.parseFasta(GBJ).features.filter(x => x.type === 'CDS')[0];
  ok(cdsJ.length === 20, 'join 的 CDS 长度 = 各段之和（20 而非 49）', cdsJ.length);
}

console.log('\n--- 17b. 虚拟克隆：末端与连接 ---');
{
  const ends = n => C.enzymeEnds(E.byName(n));

  // 末端形态必须与酶学常识一致
  const eEco = ends('EcoRI');
  ok(eEco.type === '5' && eEco.overhang === 'AATT', 'EcoRI → 5′-AATT', eEco.type + ' ' + eEco.overhang);
  const ePst = ends('PstI');
  ok(ePst.type === '3' && ePst.overhang === 'TGCA', 'PstI → 3′-TGCA', ePst.type + ' ' + ePst.overhang);
  const eSma = ends('SmaI');
  ok(eSma.type === 'blunt' && eSma.length === 0, 'SmaI → 平端');
  const eNot = ends('NotI');
  ok(eNot.type === '5' && eNot.overhang === 'GGCC', 'NotI → 5′-GGCC', eNot.overhang);

  // 兼容性：同尾酶与自身
  ok(C.endsCompatible(ends('EcoRI'), ends('EcoRI')) === true, 'EcoRI 自连可行');
  ok(C.endsCompatible(ends('BamHI'), ends('BglII')) === true, '同尾酶 BamHI × BglII 可连');
  ok(C.endsCompatible(ends('SalI'), ends('XhoI')) === true, '同尾酶 SalI × XhoI 可连');
  ok(C.endsCompatible(ends('SmaI'), ends('EcoRV')) === true, '不同平端酶可连');
  ok(C.endsCompatible(ends('EcoRI'), ends('BamHI')) === false, 'EcoRI × BamHI 不可连');
  ok(C.endsCompatible(ends('PstI'), ends('EcoRI')) === false, '3′ 突出 × 5′ 突出 不可连');
  ok(C.endsCompatible(null, ends('EcoRI')) === false, '缺末端时不误判为可连');

  // 端到端：载体切成骨架 + 插入片段，连接后得到重组序列
  const H = 'GAATTC';
  const vec = 'A'.repeat(500) + H + 'C'.repeat(400) + H + 'G'.repeat(300);   // 环状，2 个 EcoRI
  const ins = H + 'T'.repeat(500) + H;                                      // 线性，2 个 EcoRI
  const eco = E.byName('EcoRI');
  const vFrags = C.digestOneDetailed(vec, eco, true);
  const iFrags = C.digestOneDetailed(ins, eco, false);

  ok(vFrags.length === 2, '环状载体切出 2 个片段', vFrags.map(f => f.length).join('+'));
  ok(vFrags.reduce((a, f) => a + f.length, 0) === vec.length, '载体片段合计 = 原长');
  ok(vFrags.every(f => f.leftEnd && f.rightEnd), '环状片段两端都带末端信息');
  ok(iFrags.length === 3, '线性插入切出 3 个片段', iFrags.map(f => f.length).join('+'));
  ok(iFrags[0].leftEnd === null, '线性首片段左端为空（无配对对象）');
  ok(iFrags[2].rightEnd === null, '线性末片段右端为空');

  const backbone = vFrags.reduce((a, b) => (a.length > b.length ? a : b));
  const insert = iFrags.slice().sort((a, b) => b.length - a.length)[0];
  const lig = C.ligateInsert(backbone, insert);
  ok(lig.ok === true, '同酶切出的骨架与插入可连接', lig.reason);
  ok(lig.size === backbone.length + insert.length, '重组大小 = 骨架 + 插入',
    lig.size + ' = ' + backbone.length + ' + ' + insert.length);
  ok(lig.seq.slice(0, backbone.length) === backbone.seq, '重组序列前段 = 骨架');
  ok(lig.seq.slice(backbone.length) === insert.seq, '重组序列后段 = 插入');
  // 连接点应恢复出完整的 EcoRI 识别序列
  const j1 = lig.seq.slice(backbone.length - 1, backbone.length + 5);
  ok(j1 === 'GAATTC', '连接点 1 恢复出 EcoRI 位点', j1);
  const j2 = lig.seq.slice(lig.seq.length - 1) + lig.seq.slice(0, 5);
  ok(j2 === 'GAATTC', '连接点 2（跨接缝）恢复出 EcoRI 位点', j2);

  // 不兼容组合必须被拒绝
  const bad = C.ligateInsert(backbone,
    { seq: 'B'.repeat(100).replace(/B/g, 'A') + 'A'.repeat(0), leftEnd: ends('BamHI'), rightEnd: ends('BamHI') });
  ok(bad.ok === false, '末端不匹配时拒绝连接并给出原因', bad.reason.slice(0, 20));
  ok(C.ligateInsert(null, insert).ok === false, '缺骨架时安全返回');
  ok(C.ligateInsert(backbone, null).ok === false, '缺插入时安全返回');
  // 简并突出端不应被误判为可连
  const ava = ends('AvaI');
  ok(ava.degenerate === true, '简并酶的突出端被标记为无法判定', ava.overhang);
  ok(C.endsCompatible(ava, ava) === false, '简并突出端不误判为可连');
}

/* ---------- 18. 蛋白理化性质 ---------- */
console.log('\n== 18. 蛋白理化性质 ==');
{
  // 手工锚点：PP 二肽 DIWV[P][P]=20.26 → 10/2×20.26 = 101.30
  const pp = C.proteinProperties('PP');
  ok(Math.abs(pp.instability - 101.3) < 0.01, 'PP 二肽不稳定指数 = 101.30（手工锚点）', pp.instability);
  ok(pp.unstable === true, 'PP 被判为不稳定');
  ok(C.proteinProperties('KKKKK').pI > 10, 'poly-K 的 pI > 10（强碱性）');
  ok(C.proteinProperties('DDDDD').pI < 4, 'poly-D 的 pI < 4（强酸性）');
  // GRAVY 手算：'AI' = (1.8+4.5)/2 = 3.15
  ok(Math.abs(C.proteinProperties('AI').gravy - 3.15) < 1e-9, 'AI 的 GRAVY = 3.15（手算锚点）');
  // 脂肪族指数手算：'AAVV' = 50 + 2.9×50 = 195
  ok(Math.abs(C.proteinProperties('AAVV').aliphatic - 195) < 1e-9, 'AAVV 脂肪族指数 = 195（手算锚点）');
  // 消光系数：'WY' = 5500+1490 = 6990
  const wy = C.proteinProperties('WY');
  ok(wy.extinctionOx === 6990, 'WY 消光系数（还原）= 6990', wy.extinctionOx);
  ok(C.proteinProperties('C').extinctionOx - C.proteinProperties('C').extinctionRed === 125, '单个 C 贡献 125（氧化-还原差）');
  // 净电荷符号：poly-K 在 pH7 为正，poly-D 为负
  ok(C.proteinProperties('KKKKK').charge7 > 0, 'poly-K 在 pH 7 带正电');
  ok(C.proteinProperties('DDDDD').charge7 < 0, 'poly-D 在 pH 7 带负电');
  // 边界
  ok(C.proteinProperties('') === null, '空序列返回 null');
  ok(C.proteinProperties('XZBJ123') === null, '全非常规字母返回 null');
  ok(C.proteinProperties('MXZJ') !== null && C.proteinProperties('MXZJ').length === 1, '混入 X/Z/J 时只保留标准残基');
  ok(C.proteinProperties('A') !== null && C.proteinProperties('A').length === 1, '单残基可分析');
  ok(C.proteinProperties('MKTA*') !== null && C.proteinProperties('MKTA*').length === 4, '终止符被忽略');
  // pI 一定落在 0–14
  const rnd = C.proteinProperties('MKTAYIAKQRQISFVKSHFSRQLEDLRQFIERTKKLD');
  ok(rnd.pI >= 0 && rnd.pI <= 14, 'pI 落在 0–14 之间', rnd.pI);
  ok(rnd.mw > 0 && rnd.gravy !== undefined, '长蛋白各字段齐全');
}

/* ---------- 19. 反向翻译 / 密码子优化 ---------- */
console.log('\n== 19. 反向翻译 ==');
{
  const rt = C.reverseTranslate('MKTAYIAKQRQISFVKSHFSRQ');
  ok(rt.backTranslateOk === true, '回译与原蛋白完全一致');
  ok(rt.dna.length === 66, '22 aa → 66 nt', rt.dna.length);
  ok(rt.dna.slice(0, 3) === 'ATG', 'M 固定为 ATG');
  ok(C.translate(rt.dna) === 'MKTAYIAKQRQISFVKSHFSRQ', 'translate(结果) = 原序列');
  const rtStop = C.reverseTranslate('MK*');
  ok(rtStop.dna.slice(-3) === 'TAA', '终止密码子用 TAA', rtStop.dna);
  const rtBad = C.reverseTranslate('MKZXJ');
  ok(rtBad.skipped === 3, '非法残基被跳过并计数', rtBad.skipped);
  ok(rtBad.dna === 'ATGAAA', '跳过非法残基后只翻译合法部分', rtBad.dna);
  ok(C.reverseTranslate('').dna === '', '空输入返回空 DNA');
}

/* ---------- 20. 连接反应用量 ---------- */
console.log('\n== 20. 连接用量 ==');
{
  const lig = C.ligationAmounts(5000, 1000, 50, 3);
  ok(Math.abs(lig.insertNg - 30) < 1e-9, '5 kb 载体 50 ng + 1 kb 插入 3:1 → 30 ng（手算锚点）', lig.insertNg);
  ok(Math.abs(lig.insertPmol - lig.vectorPmol * 3) < 1e-12, '摩尔比严格成立');
  ok(C.ligationAmounts(0, 1000, 50, 3) === null, '载体长度为 0 返回 null');
  ok(C.ligationAmounts(5000, 1000, -5, 3) === null, '负用量返回 null');
  ok(C.ligationAmounts(5000, 1000, 50, 0) === null, '摩尔比为 0 返回 null');
  const lig2 = C.ligationAmounts(3000, 3000, 100, 1);
  ok(Math.abs(lig2.insertNg - 100) < 1e-9, '等长 1:1 → 等质量', lig2.insertNg);
}

/* ---------- 21. 定点突变引物 ---------- */
console.log('\n== 21. 定点突变引物 ==');
{
  const tpl = C.randomSeq(3000, 50);
  const r = C.designMutPrimers(tpl, 1500, 'G', { circular: true });
  ok(r.ok === true, '点突变设计成功');
  ok(r.fwd.charAt(0) === 'G', '正向引物 5′ 端即突变碱基');
  ok((tpl + tpl).indexOf(r.fwdAnneal) >= 0, '正向退火区与模板下游完全一致');
  ok((tpl + tpl).indexOf(C.revComp(r.revAnneal)) >= 0, '反向退火区反向互补于模板上游');
  ok(r.tmF >= 55 && r.tmR >= 55, '退火区 Tm 达到可用范围', r.tmF.toFixed(1) + '/' + r.tmR.toFixed(1));
  ok(r.mutDesc === tpl.charAt(1499) + '1500G', '突变描述正确', r.mutDesc);
  ok(r.productLen === 3000, '点突变产物长度不变');
  // 缺失
  const del = C.designMutPrimers(tpl, 100, '', { oldLen: 6, circular: true });
  ok(del.ok && del.productLen === 3000 - 6, '缺失 6 bp 产物长度正确', del.productLen);
  ok(del.fwd === del.fwdAnneal, '缺失突变时正向引物无 5′ 尾巴');
  // 插入
  const ins = C.designMutPrimers(tpl, 200, 'AAAAAA', { oldLen: 0, circular: true });
  ok(ins.ok && ins.productLen === 3006, '插入 6 bp 产物长度正确', ins.productLen);
  ok(ins.fwd.slice(0, 6) === 'AAAAAA', '插入序列位于正向 5′ 端');
  // 边界
  ok(C.designMutPrimers(tpl, 0, 'G', {}).ok === false, '位置 0 被拒绝');
  ok(C.designMutPrimers(tpl, 3001, 'G', {}).ok === false, '位置超出长度被拒绝');
  ok(C.designMutPrimers('', 1, 'G', {}).ok === false, '空模板被拒绝');
  ok(C.designMutPrimers(tpl, 1500, 'G', { oldLen: 4000, circular: true }).ok === false,
    '环状下被替换碱基数超过模板全长被拒绝（避免产物长度变负）');
  // 线性模板靠近末端也能设计（退火区在界内）
  const edge = C.designMutPrimers(tpl, 2, 'A', { circular: false });
  ok(edge.ok === true, '线性模板近末端点突变仍可设计');
  // 模拟验证：突变引物在突变模板上的退火区仍能配对（退火区不含突变碱基）
  const mutated = tpl.slice(0, 1499) + 'G' + tpl.slice(1500);
  ok(mutated.indexOf(r.fwdAnneal) >= 0, '退火区在突变模板上依然匹配');
}

/* ---------- 22. Motif 搜索 ---------- */
console.log('\n== 22. Motif 搜索 ==');
{
  const hits = C.findMotifs('ATGAATAAAGGAGGTATATAAT', 'AATAAA', {});
  ok(hits.length === 1 && hits[0].pos === 4 && hits[0].strand === '+', '正链单命中定位正确', JSON.stringify(hits[0]));
  const both = C.findMotifs('AGGAGGTTTTAGGAGG', 'AGGAGG', { bothStrands: true });
  ok(both.length === 2, '双链搜索命中 2 处');
  // 负链命中映射回正链坐标：序列含 CCTCCT（AGGAGG 的反向互补）
  const neg = C.findMotifs('TTTTCCTCCTTTTT', 'AGGAGG', { bothStrands: true });
  ok(neg.length === 1 && neg[0].strand === '-' && neg[0].pos === 5, '负链命中坐标映射正确', JSON.stringify(neg[0]));
  // 简并碱基：GTYRAC 中 Y=C/T、R=A/G；GTTAAC 与 GTTGAC 均匹配，GTGCAC（第3位G不满足Y）为对照
  const deg = C.findMotifs('GTTAACGTTGACGTGCAC', 'GTYRAC', {});
  ok(deg.length === 2, 'IUPAC 简并模式 GTYRAC 命中 2 处（第 3 个为对照）', deg.length);
  // 重叠命中
  const ovl = C.findMotifs('AAAAAA', 'AAA', {});
  ok(ovl.length === 4, '重叠命中全部报出', ovl.length);
  // 环状跨接缝
  const wrap = C.findMotifs('ATAAAGG', 'GGAT', { circular: true });
  ok(wrap.length === 1 && wrap[0].pos === 6, '环状跨接缝命中', JSON.stringify(wrap[0]));
  const noWrap = C.findMotifs('ATAAAGG', 'GGAT', { circular: false });
  ok(noWrap.length === 0, '线性模式漏报跨接缝（对照）');
  ok(C.findMotifs('ACGT', '', {}).length === 0, '空模式返回空数组');
  ok(C.findMotifs('', 'AAA', {}).length === 0, '空序列返回空数组');
  // 跨接缝标记
  const wrapHit = C.findMotifs('ATAAAGG', 'GGAT', { circular: true })[0];
  ok(wrapHit.wrap === true, '正链跨接缝命中带 wrap 标记');
  ok(C.findMotifs('ATAAAGGATAAAGG', 'GGAT', { circular: false }).every(h => !h.wrap), '线性命中不带 wrap 标记');
  // 负链跨接缝：seq[4..7]+seq[1]='TAGCA'，其反向互补为 TGCTA → 应在 pos=4 报一个负链命中
  const negWrap = C.findMotifs('AGGTAGC', 'TGCTA', { circular: true, bothStrands: true });
  ok(negWrap.length === 1 && negWrap[0].strand === '-' && negWrap[0].pos === 4 && negWrap[0].wrap === true,
    '负链跨接缝命中坐标正确（早期版本会算出 0 而丢命中）', JSON.stringify(negWrap.map(h => h.pos + h.strand)));
  // 结果上限
  const bigSeq = C.randomSeq(30000, 50);
  const capped = C.findMotifs(bigSeq, 'N', { maxResults: 1000 });
  ok(capped.length === 1000 && capped.truncated === true, 'maxResults 生效并标记截断');
  ok(C.findMotifs(bigSeq, 'N', {}).truncated === false, '未触上限时不标记截断');
  ok(Object.keys(capped[0]).indexOf('truncated') < 0, 'truncated 不污染结果对象');
}

/* ---------- 23. 双酶切兼容性 ---------- */
console.log('\n== 23. 双酶切兼容性 ==');
{
  const bam = E.byName('BamHI'), bgl = E.byName('BglII'), eco = E.byName('EcoRI'), hind = E.byName('HindIII');
  const seq = 'AAAAAAGGATCCAAAAAAAAAAAGATCTAAAAA';
  const r = C.doubleDigestSummary(seq, bam, bgl, false);
  ok(r.cutsA === 1 && r.cutsB === 1, '各 1 个切点');
  ok(r.compatible === true, 'BamHI×BglII 同尾互补被识别');
  ok(r.fragments.length === 3, '双酶切产生 3 个片段', JSON.stringify(r.fragments));
  const r2 = C.doubleDigestSummary(seq, eco, hind, false);
  ok(r2.cutsA === 0 && r2.cutsB === 0, '无切点酶的片段报告', JSON.stringify(r2.fragments));
  ok(r2.compatible === false, 'EcoRI×HindIII 不同尾');
  const same = C.doubleDigestSummary(seq, bam, bam, false);
  ok(same.notes.some(function (n) { return n.indexOf('相同') >= 0; }), '同酶双选给出提示');
  ok(C.doubleDigestSummary(seq, bam, null, false) === null, '缺酶返回 null');
  // 联合片段与 digest() 一致
  const dBoth = C.digest(seq, [bam, bgl], false);
  ok(JSON.stringify(r.fragments) === JSON.stringify(dBoth.fragments), '联合片段与 digest() 完全一致');
}

console.log('\n' + (fails === 0 ? '全部通过 ✔' : fails + ' 项失败 ✘'));
process.exit(fails === 0 ? 0 : 1);
