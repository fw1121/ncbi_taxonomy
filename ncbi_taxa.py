#!/bin/env python 
import sys
import os
sys.path.insert(0, os.path.expanduser("~/.pythonlib/"))

from collections import defaultdict, deque
from itertools import permutations
from argparse import ArgumentParser
from string import strip
import logging as log

import operator
import sqlite3

from ete2a1 import PhyloTree

__DESCRIPTION__ = """ 
Query ncbi taxonomy using a local DB

(jhcepas [at] gmail.com) 08 Aug 2011
"""

log.basicConfig(level=log.INFO, \
                    format="%(levelname)s - %(message)s" )

module_path = os.path.split(__file__)[0]
c = sqlite3.connect(os.path.join(module_path, 'taxa.sqlite'))


def get_sp_lineage(taxid):
    result = c.execute('select track from species where taxid=%s' %taxid)
    raw_track = result.fetchone()
    if not raw_track:
        raise ValueError("%s taxid not found" %taxid)
    track = map(int, raw_track[0].split(","))
    return reversed(track)
    
def get_taxid_translator(taxids):
    def get_name(taxid):
        result = c.execute('select spname from species where taxid=%s' %sp)
        try:
            return result.fetchone()[0]
        except TypeError:
            raise ValueError("%s taxid not found" %taxid)
    id2name = {}
    for sp in taxids: 
        id2name[sp] = get_name(sp)
    return id2name

def get_name_translator(names):
    name2id = {}
    for sp in names: 
        result = c.execute('select taxid from species where spname="%s"' %sp)
        try:
            name2id[name] = result.fetchone()[0]
        except TypeError:
            raise ValueError("[%s] name not found" %sp)
    return name2id

def translate_to_names(taxids):
    def get_name(taxid):
        result = c.execute('select spname from species where taxid=%s' %taxid)
        try:
            return result.fetchone()[0]
        except TypeError:
            raise ValueError("%s taxid not found" %taxid)
    id2name = {}
    names = []
    for sp in taxids: 
        names.append(id2name.setdefault(sp, get_name(sp)))
    return names

def get_topology(taxids):
    sp2track = {}
    elem2node = {}
    for sp in taxids:
        track = deque()
        for elem in get_sp_lineage(sp):
            node = elem2node.setdefault(elem, PhyloTree())
            node.name = str(elem)
            track.append(node)
        sp2track[sp] = track

    # generate parent child relationships
    for sp, track in sp2track.iteritems():
        parent = None
        for elem in track:
            if parent and elem not in parent.children:
                parent.add_child(elem)
            parent = elem
    root = elem2node[1]

    #remove onechild-nodes
    clean=1
    if clean:
        for n in root.get_descendants():
            if len(n.children) == 1: 
                n.delete()
    if len(root.children) == 1:
        return root.children[0].detach()
    else:
        return root

def annotate_raw_tree(t):
    tax2name = get_taxid_translator([n.name for n in t.traverse()])
    leaf2track = dict([ [n.name, get_sp_lineage(n.name)] for n in t.get_leaves()])
    for n in t.traverse():
        n.add_feature("taxid", n.name)
        n.name = tax2name[n.taxid]
        if n.is_leaf():
            n.lineage = leaf2track[n.taxid]
            n.named_lineage = translate_to_names(leaf2track[n.taxid])
    return t

def test():
    # TESTS
    get_sp_lineage("9606")
    t = get_topology([9913,31033,7955,9606,7719,9615,44689,10116,7227,9031,13616,7165,8364,99883,10090,9598])
    annotate_raw_tree(t)
    print t.get_ascii(show_internal=True, compact=False)
    t.show()


if __name__ == "__main__":
    parser = ArgumentParser(description=__DESCRIPTION__)
    # name or flags - Either a name or a list of option strings, e.g. foo or -f, --foo.
    # action - The basic type of action to be taken when this argument is encountered at the command line. (store, store_const, store_true, store_false, append, append_const, version)
    # nargs - The number of command-line arguments that should be consumed. (N, ? (one or default), * (all 1 or more), + (more than 1) )
    # const - A constant value required by some action and nargs selections. 
    # default - The value produced if the argument is absent from the command line.
    # type - The type to which the command-line argument should be converted.
    # choices - A container of the allowable values for the argument.
    # required - Whether or not the command-line option may be omitted (optionals only).
    # help - A brief description of what the argument does.
    # metavar - A name for the argument in usage messages.
    # dest - The name of the attribute to be added to the object returned by parse_args().

    parser.add_argument("-t", "--taxid", dest="taxid", nargs="+",  
                        type=int, 
                        help="""taxids (space sparated)""")

    parser.add_argument("-tf", "--taxid_file", dest="taxid_file",   
                        type=str, 
                        help="""file containing a list of taxids (one per line)""")


    parser.add_argument("-n", "--name", dest="names", nargs="+",  
                        type=str, 
                        help="""species or taxa names (comma separated)""")

    parser.add_argument("-nf", "--names_file", dest="names_file",   
                        type=str, 
                        help="""file containing a list of taxids (one per line)""")


    args = parser.parse_args()
    all_taxids = []
    if args.taxid_file:
        all_taxids.extend(map(strip, open(args.taxid_file, "rU").read().split("\n")))
    if args.taxid:
        all_taxids.extend(args.taxid)
    if all_taxids:
        log.info("Dumping %d taxid translations:" %len(all_taxids))
        all_taxids = set(all_taxids)
        all_taxids.discard("")
        for taxid, name in get_taxid_translator(all_taxids).iteritems():
            lineage = get_sp_lineage(taxid)
            named_lineage = ','.join(translate_to_names(lineage))
            lineage = ','.join(map(str, lineage))
            print "\t".join(map(str, [taxid, name, named_lineage, lineage ]))
         
    all_names = []
    if args.names_file:
        all_names.extend(map(strip, open(args.names_file, "rU").read().split("\n")))
    if args.names:
        all_names.extend(map(strip, " ".join(args.names).split(",")))
    if all_names:
        log.info("Dumping name translations:")
        for name, taxid in get_name_translator(all_names).iteritems():
            print "\t".join(map(str, [name, taxid]))
