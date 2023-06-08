from darpins import *
from subprocess import Popen
from argparse import ArgumentParser
import levenshtein as Clevenshtein
import datetime
import random
import pprint

combined_seqs = ''
ncombined = 0

def unlist_alphabet(ab):
    li = []
    for g in ab:
        li.extend(ab[g])
    return li
        
def validate_mask_lengths(seq,masks):
    for mask in masks:
        if len(masks[mask]) != len(seq):
            sys.stderr.write("CRITICAL: mask '%s' does not match sequence length\n" % mask)
            sys.exit(1)
            
def validate_mask_names(names,masks):
    for name in names:
        if name not in masks:
            sys.stderr.write("CRITICAL: name '%s' does not match any mask\n" % name)
            sys.exit(1)

def validate_alphabet_scheme(scheme):
    if scheme and scheme not in ALPHABETS:
        sys.stderr.write("CRITICAL: alphabet scheme '%s' does not exist\n" % scheme)
        sys.exit(1)

def parse_args():
    parser = ArgumentParser()

    parser.add_argument("-p", "--parent", dest="parent", required=False,
                        help="Use sequence from database from uid")
    parser.add_argument("-n", "--ndesign", dest="ndesign", type=int,
                        help="Will create <n> designs from template")
    parser.add_argument("--ntarget", dest="ntarget", type=int, required=False, default=0,
                        help="Will stop iterating after DB hits <ntarget> designs")

    parser.add_argument("-r", "--nres", dest="nres", default=0, type=int,
                        help="Will mutate at most <n> residues")
    parser.add_argument("--permutate", dest="permutate", action="store_true", default=False,
                        help="Permutates residues instead of mutating them")
    parser.add_argument("--force_seq", dest="force_seq", default="",
                        help="Forces mutation to the provided sequence")
    parser.add_argument("--force_parent", dest="force_parent", default="",
                        help="Forces creation of new parent from the provided sequence")
    parser.add_argument("--alphabet_scheme", dest="alphabet_scheme", default='',
                        help="Uses a reduced alphabet when mutating residues")
    parser.add_argument("--alphabet_change", dest="alphabet_change", action="store_true",
                        help="Forces a change in the alphabet group when mutating")
    
    parser.add_argument("-m", "--mapfile", dest="mapfile",
                      help="Masking file to map residues to mutate")
    parser.add_argument("-v", "--variable", nargs='+', dest="variable", default=['variable'],
                        help="Apply variable masking positions defined from mapfile")
    parser.add_argument("-f", "--fixed", nargs='+', dest="fixed", default=['back_res','below_res'],
                        help="Apply fixed masking positions defined from mapfile")
    parser.add_argument("--test", dest="test", default=False, action="store_true",
                        help="Use test collections from MongoDB")
    parser.add_argument("--no_insert", dest="insert", default=True, action="store_false",
                        help="Do not add entries to the database")
    
    parser.add_argument("-l", "--lev_threshold", dest="lev_threshold", default=0, type=int,
                      help="Adjust tolerance for sequence divergence")
    parser.add_argument("--lev_python", dest="lev_python", default=False, action="store_true",
                      help="Use python version of the levenshtein implementation (slower)")
    parser.add_argument("--div_python", dest="div_python", default=False, action="store_true",
                      help="Use python version to calculate pairwise levenshtein (slower)")
    parser.add_argument("--no_filter_type", dest="filter_type", default=True, action="store_false",
                      help="Do not filter by type when designing")
    parser.add_argument("--no_filter_parent", dest="filter_parent", default=True, action="store_false",
                      help="Do not filter by parent when designing")
    parser.add_argument("-s", "--seed", dest="seed", default=None,
                      help="Sets the seed of the design")

    return parser.parse_args()

def apply_mask(li,mask,val=1):
    for i in range(0,len(li)):
        if mask[i]:
            li[i] = val

def apply_masks(li,masks,names,val=1):
    for name in names:
        apply_mask(li,masks[name],val=val)

def apply_sequence_mask(seq,mask,fun):
    for i in range(0,len(seq)):
        if mask[i]:
            if fun == 'upper':
                seq[i] = seq[i].upper()
            elif fun == 'lower':
                seq[i] = seq[i].lower()
    
def mutate_residue(from_res,args,to=''):
    if to: return to
    if args.alphabet_scheme:
        ab = dict(ALPHABETS[args.alphabet_scheme])
        g0 = REVERSED_ALPHABETS[args.alphabet_scheme][from_res]
        if args.alphabet_change:
            del ab[g0]
        li = unlist_alphabet(ab)
    else:
        li = list(MUTATABLE_AA)
    # do not mutate to itself
    if from_res in li:
        li.remove(from_res)
    ipos = random.randint(0,len(li)-1)
    return li[ipos]

def levenshtein(s, t):
    if s == "":
        return len(t)
    if t == "":
        return len(s)
    cost = 0 if s[-1] == t[-1] else 1
    
    i1 = (s[:-1], t)
    if not i1 in memo:
        memo[i1] = levenshtein(*i1)
    i2 = (s, t[:-1])
    if not i2 in memo:
        memo[i2] = levenshtein(*i2)
    i3 = (s[:-1], t[:-1])
    if not i3 in memo:
        memo[i3] = levenshtein(*i3)
    res = min([memo[i1]+1, memo[i2]+1, memo[i3]+cost])
    
    return res
    
def mutate_towards_alphabet(seq,seq2,rab):
    for i in range(0,len(seq)):
        if seq2[i].isupper() and \
           rab[seq[i]] == rab[seq2[i]]:
            seq[i] = seq2[i]
    
def mutate_towards_alphabet_first(seq,rab,ab):
    for i in range(0,len(seq)):
        if seq[i].isupper():
            seq[i] = ab[rab[seq[i]]][0]
    
def is_divergent(seq,masks,designs,args):
    if not args.div_python:
        return is_divergent_full(seq,masks,designs,args)
    # print("calculating divergence...")
    min_l = 100000
    # max_l = 0
    seq_ = [ seq[i] for i in range(len(seq)) if masks['merged'][i] ]
    filterdata = {}
    if args.filter_type:
        filterdata.update({ 'type': 'permutate' if args.permutate else 'random' })
    if args.filter_parent:
        filterdata.update({ 'parent' : args.parent })
    lev_fun = levenshtein if args.lev_python else Clevenshtein.levenshtein
    for d in designs.find(filterdata):
        # provide a mutated sequence with same aa from group
        # to consider to take into account convergence based on alphabet
        global memo
        memo = {}
        seq2_ = [ d['seq'][i] for i in range(len(d['seq'])) if masks['merged'][i] ]
        if args.alphabet_scheme:
            mutate_towards_alphabet(seq2_,seq_,REVERSED_ALPHABETS[args.alphabet_scheme])
        l = lev_fun(''.join(seq_),''.join(seq2_))
        # if l > max_l: max_l = l
        if l < min_l:
            min_l = l
            if min_l < args.lev_threshold: break
    # print("min_levenshtein", min_l)
    return min_l >= args.lev_threshold
    
def is_divergent_full(seq,masks,designs,args):
    global combined_seqs;
    global ncombined;
    seq_ = [ seq[i] for i in range(len(seq)) if masks['merged'][i] ]
    if args.alphabet_scheme:
        mutate_towards_alphabet_first(seq_,REVERSED_ALPHABETS[args.alphabet_scheme],
                                      ALPHABETS[args.alphabet_scheme])
    filterdata = { 'type' : 'permutate' if args.permutate else 'random',
                   'parent' : args.parent
    }
    count = designs.count(filterdata)
    if args.ntarget and count >= args.ntarget:
        sys.exit(1)
    elif count > ncombined:
        for d in designs.find(filterdata)[ncombined:]:
            lseq_ = [ d['seq'][i] for i in range(len(d['seq'])) if masks['merged'][i] ]
            if args.alphabet_scheme:
                mutate_towards_alphabet_first(lseq_,REVERSED_ALPHABETS[args.alphabet_scheme],
                                              ALPHABETS[args.alphabet_scheme])
            combined_seqs += ''.join(lseq_)
            ncombined += 1
    return Clevenshtein.multiple_ndiff(''.join(seq_),combined_seqs,len(seq_),args.lev_threshold)
    
def mutate_sequence(seq,mutmask,indexes,args,force_seq=[]):
    # find position to mutate
    mutpos = random.randint(0,len(indexes)-1)
    i = indexes[mutpos]
    to = force_seq[i] if force_seq else ''
    seq[i] = mutate_residue(seq[i],args,to=to)
    mutmask[i] = 1
    return mutpos

def permutate_sequence(seq,mutmask,indexes,args,force_seq=[]):
    # find position to mutate
    permpos1 = random.randint(0,len(indexes)-1)
    while True:        
        permpos2 = random.randint(0,len(indexes)-1)
        if permpos2 != permpos1: break
    i = indexes[permpos1]
    j = indexes[permpos2]
    seqi = seq[i]
    seq[i] = seq[j]
    seq[j] = seqi
    mutmask[i] = 1
    mutmask[j] = 1
    return permpos1, permpos2
    
def design_sequence(seq,masks,args,force_seq=''):
    seq = [ r.lower() for r in seq ]
    for i in range(0,len(masks['merged'])):
        if masks['merged'][i]:
            seq[i] = seq[i].upper()
    nres = len([c for c in seq if c.isupper()])
    if args.nres > 0 and args.nres < nres:
        nres = args.nres
    mutmask = [0 for c in seq]
    force_seq = [ r.lower() for r in force_seq ]
    indexes = []
    for i in range(0,len(seq)):
        if seq[i].isupper():
            indexes.append(i)
    # print([ r for r in seq if r.isupper() ])
    while nres > 0:
        if args.permutate:
            permpos1, permpos2 = permutate_sequence(seq,mutmask,indexes,args,force_seq=force_seq)
            # print("permutated %d and %d" % (permpos1,permpos2))
        else:
            mutpos = mutate_sequence(seq,mutmask,indexes,args,force_seq=force_seq)
            # print("mutated %d" % (mutpos))
            # drop index to not mutate same residue twice
            indexes.pop(mutpos)
        nres -= 1
    # print([ r for r in seq if r.isupper() ])
    return seq, mutmask

def design_molecule(seq,masks,designs,parent,args,j,force_seq=''):
    print("generating sequence...")
    i = j
    while True:
        newseq, mutmask = design_sequence(seq,masks,args,force_seq=force_seq)
        lseq = ''.join(newseq)
        useq = lseq.upper()
        uid = hash(useq)
        if not designs.find_one({ 'uid': uid }) and \
           (force_seq or is_divergent(newseq,masks,designs,args)): # does it respect Levenshtein distances
            break
        i += 1
        if i % 1000 == 0:
            print("%d designs scanned" % i)
        # the next line is executed when a force_seq does not go through because it already exists in the DB
        if force_seq: return
            
    if force_seq:
        type = 'manual'
    elif args.permutate:
        type = 'permutate'
    else:
        type = 'random'
        
    shortuid = uid[0:N_SHORTUID_CHARS]
    d = { 'uid': uid, 'shortuid': uid[0:N_SHORTUID_CHARS], 'seq': useq,
          'mutmask': mutmask, 'varmask': masks['merged'],
          'parent': parent, 'created': datetime.datetime.utcnow(),
          'type': type,
          'alphabet_scheme': args.alphabet_scheme,
          'alphabet_change': args.alphabet_change,
          'lev_threshold': args.lev_threshold,
          'nres': args.nres
    }
    print(d)
    if args.insert: designs.insert_one(d)
    return i

def design(seq,masks,designs,parent,args,force_seq=''):
    if force_seq:
        design_molecule(seq,masks,designs,parent,args,0,force_seq=force_seq)
    else:
        j = 0
        if args.ndesign:
            for i in range(0,args.ndesign):
                j = design_molecule(seq,masks,designs,parent,args,j)


args = parse_args()
print(args)


designs = get_mongo_designs(test=args.test)

if args.force_parent:
    if not designs.find_one({ 'uid': hash_seq(args.force_parent) }):
        d = { 'uid': hash_seq(args.force_parent), 'shortuid': hash_seq(args.force_parent)[0:N_SHORTUID_CHARS], 'seq': args.force_parent, 
              'mutmask': [0 for c in args.force_parent], 'varmask': [0 for c in args.force_parent], 
              'parent': None, 'created': datetime.datetime.utcnow(), 'type': 'template' }
        if args.insert:
            designs.insert_one(d)
            designs.create_index([('uid', pymongo.ASCENDING)],unique=True)
            designs.create_index([('shortuid', pymongo.ASCENDING)],unique=True)
            designs.create_index([('type', pymongo.ASCENDING)])
            designs.create_index([('parent', pymongo.ASCENDING)])
            print("Successfully created the parent: %s" % hash_seq(args.force_parent)[0:N_SHORTUID_CHARS])
    else:
        print("Parent already exists: %s" % hash_seq(args.force_parent)[0:N_SHORTUID_CHARS])
    sys.exit(1)

print("number of designs in %scollection=" % ('test ' if args.test else ''), designs.count())

MUTATABLE_AA = unlist_alphabet(ALPHABETS['default'])


if args.force_seq or args.permutate:
    args.alphabet_scheme = ''
    args.alphabet_change = False
    if args.force_seq:
        args.lev_threshold = 0
        args.permutate = False
    
if not args.parent:
    print("ERROR: need to define a parent sequence from uid (-p)")
    sys.exit(1)



    
random.seed(args.seed)

args.parent = get_parent_uid(designs,args.parent)
if args.parent is None:
    sys.stderr.write("CRITICAL: could not find parent entry\n")
    sys.exit(1)

d = designs.find_one({'uid': args.parent })
seq = d['seq']
uid = hash_seq(seq)


masks = read_masks_from_file(args.mapfile)

validate_mask_lengths(seq,masks)
validate_mask_names(args.variable,masks)
validate_mask_names(args.fixed,masks)
masks['merged'] = [0 for c in seq]
apply_masks(masks['merged'],masks,args.variable,val=1)
apply_masks(masks['merged'],masks,args.fixed,val=0)


validate_alphabet_scheme(args.alphabet_scheme)


memo = {}
design(seq,masks,designs,uid,args,force_seq=args.force_seq)

print("TERMINATED")
