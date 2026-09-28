/* =============================================================
 * MolBench · enzymes.js
 * 常用限制性内切酶库（REBASE 命名，^ 表示切割位）
 * site 中 "^" 左侧碱基数 = 在识别序列内的切割偏移
 * ============================================================= */
(function (global) {
  'use strict';

  // [名称, 识别序列(含 ^), 常用缓冲液/供应商备注]
  var RAW = [
    ['AatII',   'GACGT^C',    'NEB'],
    ['Acc65I',  'G^GTACC',    'NEB'],
    ['AflII',   'C^TTAAG',    'NEB'],
    ['AgeI',    'A^CCGGT',    'NEB'],
    ['AluI',    'AG^CT',      'NEB'],
    ['ApaI',    'GGGCC^C',    'NEB'],
    ['AvaI',    'C^YCGRG',    'NEB · 简并'],
    ['AvaII',   'G^GWCC',     'NEB · 简并'],
    ['BamHI',   'G^GATCC',    'NEB'],
    ['BglII',   'A^GATCT',    'NEB'],
    ['BstBI',   'TT^CGAA',    'NEB'],
    ['ClaI',    'AT^CGAT',    'NEB'],
    ['DpnI',    'GA^TC',      'NEB · 需甲基化'],
    ['DraI',    'TTT^AAA',    'NEB'],
    ['EcoRI',   'G^AATTC',    'NEB'],
    ['EcoRV',   'GAT^ATC',    'NEB'],
    ['HaeIII',  'GG^CC',      'NEB'],
    ['HhaI',    'GC^GC',      'NEB'],
    ['HincII',  'GTY^RAC',    'NEB · 简并'],
    ['HindIII', 'A^AGCTT',    'NEB'],
    ['HpaI',    'GTT^AAC',    'NEB'],
    ['HpaII',   'C^CGG',      'NEB · 受甲基化影响'],
    ['KpnI',    'GGTAC^C',    'NEB'],
    ['MboI',    '^GATC',      'NEB · 不受甲基化影响'],
    ['MluI',    'A^CGCGT',    'NEB'],
    ['MspI',    'C^CGG',      'NEB · 不受甲基化影响'],
    ['NcoI',    'C^CATGG',    'NEB'],
    ['NdeI',    'CA^TATG',    'NEB'],
    ['NheI',    'G^CTAGC',    'NEB'],
    ['NotI',    'GC^GGCCGC',  'NEB'],
    ['NruI',    'TCG^CGA',    'NEB'],
    ['NsiI',    'ATGCA^T',    'NEB'],
    ['PmeI',    'GTTT^AAAC',  'NEB'],
    ['PstI',    'CTGCA^G',    'NEB'],
    ['PvuI',    'CGAT^CG',    'NEB'],
    ['PvuII',   'CAG^CTG',    'NEB'],
    ['RsaI',    'GT^AC',      'NEB'],
    ['SacI',    'GAGCT^C',    'NEB'],
    ['SacII',   'CCGC^GG',    'NEB'],
    ['SalI',    'G^TCGAC',    'NEB'],
    ['Sau3AI',  '^GATC',      'NEB · 同 MboI'],
    ['ScaI',    'AGT^ACT',    'NEB'],
    ['SmaI',    'CCC^GGG',    'NEB'],
    ['SnaBI',   'TAC^GTA',    'NEB'],
    ['SpeI',    'A^CTAGT',    'NEB'],
    ['SphI',    'GCATG^C',    'NEB'],
    ['StuI',    'AGG^CCT',    'NEB'],
    ['TaqI',    'T^CGA',      'NEB'],
    ['XbaI',    'T^CTAGA',    'NEB'],
    ['XhoI',    'C^TCGAG',    'NEB'],
    ['XmaI',    'C^CCGGG',    'NEB · SmaI 同裂酶']
  ];

  var ENZYMES = RAW.map(function (r) {
    return {
      name: r[0],
      site: r[1],
      note: r[2] || '',
      siteClean: r[1].replace('^', ''),
      cutOffset: r[1].indexOf('^'),
      length: r[1].replace('^', '').length
    };
  });

  var BY_NAME = {};
  ENZYMES.forEach(function (e) { BY_NAME[e.name] = e; });

  global.MolEnzymes = {
    list: ENZYMES,
    byName: function (n) { return BY_NAME[n]; },
    names: ENZYMES.map(function (e) { return e.name; })
  };
})(typeof window !== 'undefined' ? window : this);
