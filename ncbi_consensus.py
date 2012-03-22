#!/usr/bin/env python 
import sys
import os
from string import strip
from collections import defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter

sys.path.insert(0, "/home/jhuerta/_Devel/ete/gui")
from ete_dev import PhyloTree, Tree, SVG_COLORS, faces, treeview, NodeStyle, TreeStyle

import ncbi_query as ncbi

__DESCRIPTION__ = ("Calculates the consensus of a tree with the NCBI taxonomy."
                   " The analysis can be visualized over the tree, in"
                   " which broken clades are shown.")

name2color = {}
def ncbi_layout(node):
    global name2color
    if node.is_leaf():
        tax_pos = 10
        for tax,k in zip(node.lineage, node.named_lineage):
            f = faces.TextFace("%10s" %k, fsize=15)
            try:
                color = name2color[k]
            except KeyError:
                name2color[k] = color = treeview.main.random_color()
         
            if hasattr(node, "broken_groups") and tax in node.broken_groups:
                f.background.color = color
                faces.add_face_to_node(f, node, tax_pos, position="aligned")
                tax_pos += 1
        f = faces.AttrFace("spname", fsize=15)
        faces.add_face_to_node(f, node, 10, position="branch-right")
    else:
        if getattr(node, "broken_groups", None):
            f = faces.TextFace("\n".join(node.broken_groups), fsize=15, fgcolor="red")
            faces.add_face_to_node(f, node, 1, position="branch-right")
            
def analyze_tracks(t):
    tax2track = {}
    tax2name = {}
    for n in t.iter_leaves():
        tax2track[n.taxid] = list(n.lineage)
        for number, text in zip(n.lineage, n.named_lineage):
            tax2name[number] = text
    # reverse the dictionary
    track2tax = defaultdict(set)
    for tax, track in tax2track.iteritems():
        for name in track:
            track2tax[name].add(tax)
    taxid2node = dict([(n.taxid,n) for n in t.iter_leaves()])

    non_mono = set()
    mono = set()
    for trackid, taxids in track2tax.iteritems():
        if len(taxids) > 2:
            affected_nodes = [taxid2node[tax] for tax in taxids]
            common = t.get_common_ancestor(affected_nodes)
            found_taxids = set([n.taxid for n in common.iter_leaves()])
            found_taxids.discard(set([None, ""]))
            if found_taxids - taxids:
                non_mono.add((trackid, frozenset(found_taxids - taxids)))
            else:
                mono.add((trackid, frozenset(taxids)))
    
    return mono, non_mono, tax2name

def get_speciation_trees(t):
    def is_dup(x):
        sp_subtotal = sum([len(n2species[_ch]) for _ch in x.children])
        return len(n2content[x]) > 1 and len(n2species[x]) != sp_subtotal

    def get_subtrees(n):
        subtrees = []

        if is_dup(n):
            for ch in n.get_children():
                ch.detach()
                subtrees.extend(get_subtrees(ch))
        else:
            nw = n.write()
            to_visit = []
            for _n in n.iter_leaves(is_leaf_fn=is_dup):
                if is_dup(_n):
                    to_visit.append(_n)
                    #_n.img_style["bgcolor"] = treeview.main.random_color()
            #n.show()
            for _n in to_visit:
                _n.detach()
                
            freaks = [_n for _n in n.traverse() if
                       len(_n.children) == 1 or
                       (not getattr(_n, "taxid", None) and not _n.children)]
            
            for s in freaks:
                s.delete(prevent_nondicotomic=True)
            
            if n.children:
                subtrees.append(n)
            for _n in to_visit:
                subtrees.extend(get_subtrees(_n))
        return subtrees

    n2content = t.get_node2content()
    n2species = {}
    dups = 0
    for n, content in n2content.iteritems():
        n2species[n] = set([_n.taxid for _n in content])
    for n in n2content:
        if is_dup(n):
            dups += 1
            #n.img_style["fgcolor"] = "red"
            #n.img_style["size"] = 40
    #t.show()
    all_subtrees = set(get_subtrees(t))
    return all_subtrees, dups

def analyze_subtrees(subtrees):
    ncbi_mistakes = 0
    valid_subtrees = 0
    broken_groups = set()
    correct_groups = set()
    broken_subtrees = 0
    total_rf = 0
    for subt in subtrees:
        subt_size = len(subt)
        if subt_size > 1:
            valid_subtrees += 1
            if reft:
                for _n in subt.iter_leaves():
                    #_n.spcode = _n.realname
                    try:
                        _n.spcode = _n.realname
                    except:
                        print subt
                rf, rf_max = subt.robinson_foulds(reft, attr_t1="spcode")
                total_rf += float(rf)/rf_max

            if args.tax_info:
                si, no, tax2name = analyze_tracks(subt)
                ncbi_mistakes += len(no)
                if no:
                    broken_subtrees += 1
                broken_groups.update([e[0] for e in no])
                correct_groups.update([e[0] for e in si])
                children = []
                for tip in subt.iter_leaves():
                    target = (t&tip.name)
                    children.append(target)
                    target.broken_groups = set([e[0] for e in no])
                # Annotate node
                t.get_common_ancestor(children).broken_groups = set([tax2name[e[0]] for e in no])
                
    return valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf

def annotate_tree_with_taxa(t, name2taxa_file):
    names2taxid = dict([map(strip, line.split("\t"))
                        for line in open(name2taxa_file)])
    for n in t.iter_leaves():
        n.add_features(taxid=names2taxid[n.name])
    ncbi.annotate_tree(t)
        
if __name__ == "__main__":
    parser = ArgumentParser(description=__DESCRIPTION__, 
                            formatter_class=RawDescriptionHelpFormatter)
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
    
    parser.add_argument("--show", dest="show_tree",
                        action="store_true", 
                        help="""Display tree after the analysis.""")
    parser.add_argument("--render", dest="render",
                        action="store_true", 
                        help="""Display tree after the analysis.""")
    
    parser.add_argument("-t", "--tree", dest="target_tree",  nargs="+",
                        type=str, required=True,
                        help="""Tree file in newick format""")

    parser.add_argument("-tf", dest="tree_list_file",
                        type=str, 
                        help="File with the list of tree files")
    
    parser.add_argument("--tax", dest="tax_info", type=str,
                        help="If the taxid attribute is not set in the"
                        " newick file for all leaf nodes, a tab file file"
                        " with the translation of name and taxid can be"
                        " provided with this option.")

    parser.add_argument("--ref", dest="ref_tree", type=str,
                        help="Uses ref tree to compute robinson foulds"
                        " distances of the different subtrees")
                       
    parser.add_argument("-o", dest="output", type=str,
                        help="Writes result into a file")
    
    args = parser.parse_args()

    if args.output:
        OUT = open(args.output, "w")
    else:
        OUT = sys.stdout

    print >>sys.stderr, "Dumping results into", OUT
    header = ("Tree".center(50), "Total dups", "Total subtrees", "Broken subtrees", "Broken NCBI clades", "mean RF", "absolute RF")
    print >>OUT, "#"+' '.join([h.center(15) for h in header])
    target_trees = []
    if args.tree_list_file:
        target_trees = [line.strip() for line in open(args.tree_list_file)]
    if args.target_tree:
        target_trees += args.target_tree
    
    for tfile in target_trees:
        print tfile
        t = PhyloTree(tfile)
        t.sort_descendants()    

        if args.ref_tree:
            reft = Tree(args.ref_tree, format=1)
        else:
            reft = None
        
        if args.tax_info:
            annotate_tree_with_taxa(t, args.tax_info)
        else:
            for n in t.iter_leaves():
                if "_" in n.realname:
                    spcode = n.realname.split("_")[1][:3]
                else:
                    spcode = n.realname[:3]
                n.add_features(taxid=spcode)

        if reft:
            rf = t.robinson_foulds(reft, attr_t1="realname")
        else:
            rf = 0

        # Split tree into species trees
        subtrees, dups =  get_speciation_trees(t.copy())
        valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf = analyze_subtrees(subtrees)
        
        iter_values = (os.path.basename(tfile).ljust(50), dups, valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf, rf)
        print >>OUT, ' '.join([str(h).center(15)  for h in iter_values])
       
        if args.show_tree or args.render:
            ts = TreeStyle()
            ts.force_topology = True
            ts.show_leaf_name = False
            ts.layout_fn = ncbi_layout 
            ts.mode = "r"
            if args.show_tree:
                t.show(tree_style=ts)
            else:
                t.render("img.svg", tree_style=ts, dpi=300)
    if args.output:
        OUT.close()
        
