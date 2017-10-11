from lingpy import *
from lingpy.evaluate.acd import bcubes, _get_cogs, _get_bcubed_score, random_cognates
from itertools import combinations
from lingpy.compare.partial import Partial


eastern = ['NorthMarquesan_38', 'Austral_128', 'Austral_1213', 
            'Tahitian_173', 'Sikaiana_243', 'Maori_85', 'Hawaiian_52',
            'Mangareva_239', 'Tuamotuan_246', 'Rapanui_264'] 

try:
    lex = LexStat('east-polynesian-6.bin.tsv', segments='segments')
except:
    wl = Wordlist('polynesian-6.tsv')
    D = {k: wl[k] for k in [x for x in wl if wl[x, 'doculect'] in eastern]}
    D[0] = [h for h in sorted(wl.header, key=lambda x: wl.header[x])]
    wl = Wordlist(D)
    wl.add_entries('ipa', 'form', lambda x: x)
    for k in wl:
        wl[k, 'segments'] = [k for k in wl[k, 'segments'] if k not in '_+']

    lex = LexStat(wl, segments='segments')
    lex.get_scorer(runs=10000)
    lex.output('tsv', filename='east-polynesian-6.bin', ignore=[])

try:
    lex2 = LexStat('east-polynesian-6.bin2.tsv', segments='segments')
except:
    
    wl = Wordlist('polynesian-6.tsv')
    D = {k: wl[k] for k in [x for x in wl if wl[x, 'doculect'] in eastern]}
    D[0] = [h for h in sorted(wl.header, key=lambda x: wl.header[x])]
    wl = Wordlist(D)
    for k in wl:
        wl[k, 'segments'] = ' '.join(ipa2tokens(wl[k, 'form'].replace(' ',
            '_').replace('G', 'g'),
            merge_vowels=False))
        cls = tokens2class(wl[k, 'segments'].split(' '), 'sca')
        if '0' in cls:
            print(wl[k, 'doculect'], cls, wl[k, 'segments'])

    lex2 = LexStat(D, segments='segments')
    lex2.get_scorer(runs=10000)
    lex2.output('tsv', filename='east-polynesian-6.bin2', ignore=[])
    

prfs = []
nums = []
for i in range(100):
    ref = 'random_'+str(i+1)
    random_cognates(lex, ref=ref, bias='lumper')
    nums += [max(lex.get_etymdict(ref=ref))]
    p, r, f = bcubes(lex, 'cogid', 'random_'+str(i+1), pprint=False)
    prfs += [[p, r, f]]
print('{0:10}  {1:.4f}  {2:.4f}  {3:.4f}'.format(
    'RANDOM', 
    sum([x[0] for x in prfs]) / 100,
    sum([x[1] for x in prfs]) / 100,
    sum([x[2] for x in prfs]) / 100))

lex.renumber("concept")
p, r, f = bcubes(lex, 'cogid', 'conceptid', pprint=False)
print('{0:10}  {1:.4f}  {2:.4f}  {3:.4f}'.format(
    'Lumpi', p, r, f))



for method, threshold in [('edit-dist', 0.75), ('turchin', 0), ('sca', 0.45), 
        ('lexstat', 0.6), ('infomap', 0.55)]:
    mid = method.replace('-', '')+'id'
    if method == 'lexstat':
        cm = 'upgma'
    if method == 'infomap':
        cm = 'infomap'
        method = 'lexstat'
    else:
        cm = 'upgma'
    lex.cluster(method=method, threshold=threshold, ref=mid,
            cluster_method=cm, restricted_chars="")
    lex2.cluster(method=method, threshold=threshold, ref=mid,
            cluster_method=cm, restricted_chars="")

    p, r, f = bcubes(lex, 'cogid', mid, pprint=False)
    print('{0:10}  {1:.4f}  {2:.4f}  {3:.4f}'.format(mid, p, r, f))
    p, r, f = bcubes(lex2, 'cogid', mid, pprint=False)
    print('{0:10}  {1:.4f}  {2:.4f}  {3:.4f}'.format(mid+'!', p, r, f))

part = Partial(lex, segments='segments')
part.partial_cluster(method='lexstat', ref='partialids', threshold=0.55)
part.add_cognate_ids('partialids', 'partid', idtype='loose')
p, r, f = bcubes(part, 'cogid', 'partid', pprint=False)
print('{0:10}  {1:.4f}  {2:.4f}  {3:.4f}'.format('Partial', p, r, f))



# "neck"
langs = ['Austral_1213', 'Hawaiian_52', 'Mangareva_239', 'Maori_85',
        'Rapanui_264']
words = [4564, 5816, 2079, 919, 1815]

for method in ['edit-dist', 'turchin', 'sca', 'lexstat']:
    matrix = [[0 for idx in words] for idxs in words] 
    for i, w1 in enumerate(words):
        for j, w2 in enumerate(words):
            if i < j:
                if method in ['sca', 'lexstat']:
                    distance = lex.align_pairs(w1, w2, method=method,
                            distance=True, return_distance=True, pprint=False,
                            mode='overlap')
                elif method == 'edit-dist':
                    distance = edit_dist(lex[w1, 'segments'], lex[w2, 'segments'],
                            normalized=True)
                else:
                    distance = turchin(lex[w1, 'segments'], lex[w2, 'segments'])
                matrix[i][j] = distance
                matrix[j][i] = distance
    tree = upgma(matrix, taxa=langs)
    print(Tree(tree).asciiArt())
    input(method)

        
    print(method)
    print(r"\tabular{ll|p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm}}")
    for i, line in enumerate(matrix):
        print('{0:20}'.format(langs[i].split('_')[0])+r' &\\sil '+
                '{0:10}'.format(''.join(lex[words[i], 'segments']))+
                '& '+'& '.join(['{0:.2f}'.format(x) for x in 
            line])+r'\\')
    print(r'\endtabular')
    print('')
    
input('weiter?')
lexstat_wins = 0
edit_wins = 0
sca_wins = 0
all_loose = 0
for concept in lex.rows:
    cogids = _get_cogs('cogid', concept, lambda x: x, lex)
    editids = _get_cogs('editdistid', concept, lambda x: x, lex)
    turchinids = _get_cogs('turchinid', concept, lambda x: x, lex)
    scaids = _get_cogs('scaid', concept, lambda x: x, lex)
    lexstatids = _get_cogs('infomapid', concept, lambda x: x, lex)
    if cogids == lexstatids:
        if cogids != editids and cogids != scaids:
            print(concept)
            idxs = lex.get_list(row=concept, flat=True)
            words = [''.join(lex[idx, 'segments']) for idx in idxs]
            langs = [lex[idx, 'doculect'] for idx in idxs]
            for a, b, c, d, e, f, g, h in zip(idxs, words, langs, cogids, editids,
                    turchinids, scaids, lexstatids):
                print('{0:5}  {1:15}  {2:20}  {3:3}  {4:3}  {5:3}  {6:3}  {7:3}'.format(
                    a, b, c, d, e ,f, g, h))
            print('---')
    elif cogids == editids:
        if cogids != lexstatids and cogids != scaids:
            print('!!!', concept)
            idxs = lex.get_list(row=concept, flat=True)
            words = [''.join(lex[idx, 'segments']) for idx in idxs]
            langs = [lex[idx, 'doculect'] for idx in idxs]
            for a, b, c, d, e, f, g, h in zip(idxs, words, langs, cogids, editids,
                    turchinids, scaids, lexstatids):
                print('{0:5}  {1:15}  {2:20}  {3:3}  {4:3}  {5:3}  {6:3}  {7:3}'.format(
                    a, b, c, d, e ,f, g, h))
            print('---')


    if cogids == lexstatids:
        if cogids != scaids and cogids != editids:
            lexstat_wins += 1
    else:
        if editids == cogids and cogids != scaids:
            edit_wins += 1
        elif scaids == cogids and cogids != editids:
            sca_wins += 1
print(lexstat_wins, edit_wins, sca_wins)
    

