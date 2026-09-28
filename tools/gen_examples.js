/* 生成示例序列：构造保证正确的 ORF / MCS，输出 assets/examples.js 与 examples/*.fasta */
const fs = require('fs');
const path = require('path');

// 固定种子随机，保证每次生成一致
let _seed = 20260914;
function rnd() {
  _seed = (_seed * 1103515245 + 12345) & 0x7fffffff;
  return _seed / 0x7fffffff;
}
function pick(arr) { return arr[Math.floor(rnd() * arr.length)]; }

const STOPS = new Set(['TAA', 'TAG', 'TGA']);
const CODONS = [];
const B = 'ACGT';
for (const a of B) for (const b of B) for (const c of B) {
  const cod = a + b + c;
  if (!STOPS.has(cod) && cod !== 'ATG') CODONS.push(cod);
}

function randSeq(n, gc) {
  let s = '';
  for (let i = 0; i < n; i++) {
    if (rnd() < gc) s += (rnd() < 0.5 ? 'G' : 'C');
    else s += (rnd() < 0.5 ? 'A' : 'T');
  }
  return s;
}

/** 生成无内部终止子的 CDS（长度 len 为 3 的倍数，含 ATG 起始 + 终止子） */
function makeCDS(aaLen) {
  let s = 'ATG';
  for (let i = 1; i < aaLen; i++) s += pick(CODONS);
  s += 'TAA';
  return s;
}

/* ---------- 示例 1：表达盒质粒 ---------- */
const promoter = 'TTGACAATTAATCATCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCACACA' + randSeq(60, 0.42);
const rbs = 'AGGAGGTAAAAAAAC';
const cds = makeCDS(148);
const mcs =
  'GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT' +
  'GGCTCGAGCGGCCGCCCATATGCCATGGATCGATCTGCAGGAATTCC';
const term = randSeq(90, 0.38) + 'AAAAAAGCGGCCGCAAAAAA';

function gcOf(s) {
  const g = (s.match(/[GC]/g) || []).length;
  return g / s.length;
}

// 在起始密码子上游插入同框终止子 TAA，避免上游 ATG 把目标 CDS 延长成更长的 ORF
// （真实基因组中这种情况很常见，但作为教学示例需要保证目标 CDS 边界清晰）
let plasmid = promoter + rbs + 'TAA' + cds + randSeq(70, 0.45) + mcs + term;
// 补足到 ~2200 bp，并让 GC 落在 48-52%
while (plasmid.length < 2200) plasmid += randSeq(120, 0.5);
// 追加一段 GC 富集区（模拟结构域 / 重复序列），便于 GC 窗口图有形状
plasmid = plasmid.slice(0, 1600) + randSeq(260, 0.72) + plasmid.slice(1600);
// 再补一段 AT 富集区
plasmid = plasmid + randSeq(220, 0.30);
plasmid = plasmid.toUpperCase();

/* ---------- 示例 2：基因组片段（含 GC 岛） ---------- */
let genome = randSeq(400, 0.40) + randSeq(500, 0.68) + randSeq(450, 0.41) + randSeq(300, 0.55);
const cds2 = makeCDS(96);
genome = genome.slice(0, 300) + 'TAA' + cds2 + genome.slice(300);
genome = genome.toUpperCase();

const EXAMPLES = [
  {
    id: 'plasmid',
    name: 'pMB-01 大肠杆菌表达盒（环状质粒）',
    desc: '含启动子 / RBS / 148 aa CDS / 多克隆位点 MCS / 终止子，含 GC 富集区与 AT 富集区。默认按环状分析。',
    circular: true,
    seq: plasmid
  },
  {
    id: 'genomic',
    name: 'gMB-02 基因组片段（含 GC 岛）',
    desc: '1650 bp 基因组序列，中部有明显 GC 岛，含一个 97 aa 的完整 CDS。用于 GC 窗口与 ORF 演示。',
    circular: false,
    seq: genome
  }
];

const root = path.join(__dirname, '..');
const outJs = '/* MolBench 内置示例序列（由 tools/gen_examples.js 生成，勿手改） */\n' +
  'window.MolExamples = ' + JSON.stringify(EXAMPLES, null, 2) + ';\n';
fs.writeFileSync(path.join(root, 'assets', 'examples.js'), outJs, 'utf8');

const fa = EXAMPLES.map(e =>
  '>' + e.id + ' | ' + e.name + ' | ' + e.seq.length + ' bp | circular=' + e.circular + '\n' +
  e.seq.replace(/(.{60})/g, '$1\n')
).join('\n');
fs.writeFileSync(path.join(root, 'examples', 'MolBench-demo.fasta'), fa, 'utf8');

// 自检输出
console.log('示例1 长度:', plasmid.length, 'GC%:', (gcOf(plasmid) * 100).toFixed(2));
console.log('示例2 长度:', genome.length, 'GC%:', (gcOf(genome) * 100).toFixed(2));
console.log('CDS 长度:', cds.length, '是否为3倍数:', cds.length % 3 === 0);
console.log('CDS 内部终止子数:', (cds.slice(0, -3).match(/TAA|TAG|TGA/g) || []).length);
for (const e of ['GAATTC', 'GGATCC', 'AAGCTT', 'CTCGAG', 'GTCGAC', 'CTGCAG', 'GGTACC', 'CCCGGG', 'GCGGCCGC', 'CATATG', 'CCATGG']) {
  console.log('  MCS 位点', e, ':', (plasmid.match(new RegExp(e, 'g')) || []).length, '处');
}
