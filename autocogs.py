from lingpy import *
from sys import argv
from lingpy.compare.sanity import mutual_coverage_check, mutual_coverage_subset
from collections import defaultdict
from lingpy.evaluate.acd import bcubes, diff

def segments():
    seq1, seq2, seq3, seq4, seq5 = "th o x t a", "thoxta", "apfəl", "tʰoxtɐ", "dɔːtər"

    print(seq1, "\t->\t", '\t'.join(ipa2tokens(seq1)))
    print(seq2, "  \t->\t", '\t'.join(ipa2tokens(seq2)))
    print(seq2, "  \t->\t", '\t'.join(ipa2tokens(seq2, semi_diacritics="h")))
    print(seq3, "  \t->\t", '\t'.join(ipa2tokens(seq3)))
    print(seq3, "  \t->\t", '\t'.join(ipa2tokens(seq3, semi_diacritics="f")))
    print(seq4, "  \t->\t", '\t'.join(ipa2tokens(seq4)))
    print(seq5, "  \t->\t", '\t'.join(ipa2tokens(seq5)))


def evaluate():

    wl = Wordlist('mikronesian-lexstat.tsv')

    for res in ['turchinid', 'scaid', 'lexstatid', 'infomap']:
        print('{0:10}\t{1[0]:.2f}\t{1[1]:.2f}\t{1[2]:.2f}'.format(
            res,
            bcubes(wl, 'cogid', res, pprint=False)
        ))
    

def classes():

    word = "θiɣatɛra"
    segs = ipa2tokens(word)
    
    # iterate over sound class models and write them in converted version 
    for m in ['dolgo', 'sca', 'asjp', 'art']:
        print(word, ' -> ', ''.join(tokens2class(segs, m)), '({0})'.format(m))


def errors():
    def check_sequence(seq):
        """Takes a segmented string as input and returns erroneously converted segments."""
        cls = tokens2class(seq, 'dolgo') # doesn't matter which model to take, all cover the same character range
        errors = defaultdict(int)
        for t, c in zip(seq, cls):
            if c == '0':
                errors[t] += 1
        return errors
    
    word = "θiɣatEra"
    seq = ipa2tokens(word)
    for error, count in check_sequence(seq).items():
        print("The symbol <{0}> occurs {1} times and is not recognized.".format(
            error, count))


def wordlist1():

    # load the wordlist
    wl = Wordlist('polynesian.tsv')
    
    # count number of languages, number of rows, number of concepts
    print("Wordlist has {0} languages and {1} concepts across {2} rows.".format(
        wl.width, wl.height, len(wl)))


def wordlist2():
    
    wl = Wordlist('polynesian.tsv')
    # get all indices for concept "hand", `row` refers to the concepts here, while `col` refers to languages
    eight = wl.get_dict(row='Eight', entry='value')
    for taxon in ['Emae_1030', 'RennellBellona_206', 'Tuvalu_753', 'Sikaiana_243', 'Penrhyn_235',  'Kapingamarangi_217']:
        print('{0:20}'.format(taxon), '  \t', ', '.join(eight[taxon]))


def coverage1():
    
    wl = Wordlist('polynesian.tsv')
    for i in range(210, 0, -1):
        if mutual_coverage_check(wl, i):
            print("Minimal mutual coverage is at {0} concept pairs.".format(i))
            break


def coverage2():
    wl = Wordlist('polynesian.tsv')
    count, results = mutual_coverage_subset(wl, 200)
    coverage, languages = results[0]
    print('Found {0} languages with an average mutual coverage of {1}.'.format(count, coverage))
    
    # write word list to file
    wl.output("tsv", filename="mikronesian", subset=True, rows=dict(doculect = "in "+str(languages)), 
              ignore='all', prettify=False)
    
    # load the smaller word list
    wl = Wordlist('mikronesian.tsv')
    
    # print basic characteristics
    print("The new word list has {0} languages and {1} concepts across {2} words.".format(
        wl.width, wl.height, len(wl)))



def alignments1():
    msa = Multiple(['β a r u', 'v a ŋ g u', 'v a l u', 'v a l u', 'v a r u', 'w a l u'])
    msa.prog_align()
    print(msa)


def alignments2():
    msa = Multiple(['β a r u', 'v a ŋ g u', 'v a l u', 'v a l u', 'v a r u', 'w a l u'])

    msa.lib_align()
    print(msa)
    

def alignments3():
    words = ['j a b l o k o', 'j a b ə l k o', 'j a b l k o', 'j a p k o']
    msa = Multiple(words)
    msa.prog_align()
    print(msa)
    print('There is {0} swap in the alignment.'.format('no' if not msa.swap_check(swap_penalty=-1) else 'a'))


def trigger_errors():
    wl = Wordlist('mikronesian.tsv')
    
    # add new column "segments" and replace data from column "tokens"
    wl.add_entries('segments', 'tokens', lambda x: ['A' if y == 'e' else y for y in x])
    
    lex = LexStat(wl, segments='segments', check=True)


def cognates_turchin():
    lex = LexStat('mikronesian.tsv', segments='tokens', check=True)
    
    # run the dolgopolsky (turchin) analysis, which is threshold-free
    lex.cluster(method='turchin')
    
    # show the cognate sets, stored in "turchinid" for the words for "Eight"
    eight = lex.get_dict(row='Eight') # get a dictionary with language as key for concept "eight"
    for k, v in eight.items():
        idx = v[0] # index of the word, it gives us access to all data
        print("{0:20} \t {1} \t{2}".format(lex[idx, 'doculect'], lex[idx, 'value'], lex[idx, 'turchinid']))


def cognates_sca():
    
    lex = LexStat('mikronesian.tsv', segments='tokens', check=True)
    lex.cluster(method='turchin')
    lex.cluster(method="sca", threshold=0.45)
   
    # show the cognate sets, stored in "turchinid" for the words for "Eight"
    eight = lex.get_dict(row='Eight') # get a dictionary with language as key for concept "eight"

    for k, v in eight.items():
        idx = v[0] 
        print("{0:20} \t {1} \t{2} \t {3} ".format(
            lex[idx, 'doculect'], 
            lex[idx, 'value'], 
            lex[idx, 'turchinid'], 
            lex[idx, 'scaid']))


def cognates_lexstat():
    
    lex = LexStat('mikronesian.tsv', segments='tokens', check=True)
    lex.get_scorer(runs=10000)

    lex.cluster(method='turchin')
    lex.cluster(method="sca", threshold=0.45)
    # show the cognate sets, stored in "turchinid" for the words for "Eight"
    eight = lex.get_dict(row='Eight') # get a dictionary with language as key for concept "eight"   
    for k, v in eight.items():
        idx = v[0] 
        print("{0:20} \t {1} \t{2} \t {3} ".format(
            lex[idx, 'doculect'], 
            lex[idx, 'value'], 
            lex[idx, 'turchinid'], 
            lex[idx, 'scaid']))
    lex.output('tsv', filename='mikronesian.bin')

    
def cognates_infomap():
    
    try:
        lex = LexStat('mikronesian.bin.tsv', segments='tokens')
    except:
        lex = LexStat('mikronesian.tsv', segments='tokens', check=True)
        lex.get_scorer(runs=10000)
        lex.output('tsv', filename='mikronesian.bin')

    lex.cluster(method="lexstat", threshold=0.55, ref="infomap")
    eight = lex.get_dict(row='Eight') # get a dictionary with language as key for concept "eight"   
    for k, v in eight.items():
        idx = v[0] 
        print("{0:20} \t {1} \t{2} \t {3} \t {4} \t {5}".format(
            lex[idx, 'doculect'], 
            lex[idx, 'value'], 
            lex[idx, 'turchinid'], 
            lex[idx, 'scaid'],
            lex[idx, 'lexstatid'],
            lex[idx, 'infomap']
        ))

    lex.output('tsv', filename='mikronesian-lexstat', ignore='all', prettify=False)


def alignments4():

    lex = LexStat('mikronesian.bin.tsv')
    alm = Alignments('mikronesian-lexstat.tsv', segments='tokens', 
            ref='infomap') # `ref` indicates the column with the cognate sets
    alm.align(method='progressive', scoredict=lex.cscorer)
    alm.output('tsv', filename='mikronesian-aligned', ignore='all',
            prettify=False)


def alignments5():

    alm = Alignments('mikronesian-aligned.tsv', segments='tokens',
        ref='infomap')
    for cog in ['1']:
        
        msa = alm.msa['infomap'][cog]
        for i, idx in enumerate(msa['ID']):
            print(
                '{0:20}'.format(msa['taxa'][i]),  
                '\t',
                alm[idx, 'concept'],
                '\t',
                '\t'.join(msa['alignment'][i])
            )

def tree():
    
    wl = Wordlist('mikronesian-lexstat.tsv')
    wl.calculate('tree', tree_calc='upgma')
    wl.output('tre', filename='mikronesian')
    print(wl.tree.asciiArt())


def distances():

    wl = Wordlist('mikronesian-lexstat.tsv')
    wl.calculate('dst', ref='infomap', mode='swadesh')
    wl.output('dst', filename='mikronesian')


def nexus():
    wl = Wordlist('mikronesian-lexstat.tsv')
    wl.output('paps.nex', filename='mikronesian', ref='infomap', missing='?')


def diffit():
    wl = Wordlist('mikronesian-lexstat.tsv')
    diff(wl, 'cogid', 'infomap', pprint=False)
    
if __name__ == '__main__':
    
    if 'diff' in argv:
        diffit()
        
    if 'evaluate' in argv:
        evaluate()

    if 'nexus' in argv:
        nexus()
        
    if 'distances' in argv:
        distances()
        
    if 'tree' in argv:
        tree()

    if 'segments' in argv:
        segments()

    if 'classes' in argv:
        classes()

    if 'errors' in argv:
        errors()

    if 'wordlist1' in argv:
        wordlist1()

    if 'wordlist2' in argv:
        wordlist2()

    if 'coverage1' in argv:
        coverage1()

    if 'coverage2' in argv:
        coverage2()

    if 'alignments1' in argv:
        alignments1()

    if 'alignments2' in argv:
        alignments2()

    if 'cognates-turchin' in argv:
        cognates_turchin()

    if 'cognates-sca' in argv:
        cognates_sca()

    if 'alignments3' in argv:
        alignments3()

    if 'alignments4' in argv:
        alignments4()

    if 'alignments5' in argv:
        alignments5()

    if 'trigger-errors' in argv:
        trigger_errors()

    if 'cognates-infomap' in argv:
        cognates_infomap()

    if 'cognates-lexstat' in argv:
        cognates_lexstat()

    if 'help' in argv:
        print("Usage python autocogs.py CMD")
