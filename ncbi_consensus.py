#!/usr/bin/env python 
import sys
import os
import cPickle
import numpy
import time
from string import strip
from collections import defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter

#try:
sys.path.insert(0, "/home/jhuerta/_Devel/ete/master")
from ete_dev import PhyloTree, Tree, SVG_COLORS, faces, treeview, NodeStyle, TreeStyle
#except ImportError:
#    sys.path.insert(0, "/users/tg/jhuerta/ete_tiles")
#    print "Using NFS ETE source"
#    from ete_dev import PhyloTree, Tree, SVG_COLORS, faces, treeview, NodeStyle, TreeStyle
   

import ncbi_query as ncbi

__DESCRIPTION__ = ("Calculates the consensus of a tree with the NCBI taxonomy."
                   " The analysis can be visualized over the tree, in"
                   " which broken clades are shown.")

try:
    name2color = cPickle.load(open("ncbi_colors.pkl"))
except Exception:
    name2color = {}
   

def npr_layout(node):
    if node.is_leaf():
        name = faces.AttrFace("name", fsize=12)
        faces.add_face_to_node(name, node, 0, position="branch-right")
        if hasattr(node, "sequence"):
            seq_face = faces.SeqFace(node.sequence, [])
            faces.add_face_to_node(seq_face, node, 0, position="aligned")

        
    if "alg_type" in node.features:
        faces.add_face_to_node(faces.AttrFace("alg_type", fsize=8), node, 0, position="branch-top")
        ttype=faces.AttrFace("tree_type", fsize=8, fgcolor="DarkBlue")
        faces.add_face_to_node(ttype, node, 0, position="branch-top")
        #ttype.background.color = "DarkOliveGreen"
        node.img_style["size"] = 20
        node.img_style["fgcolor"] = "red"
        
    if "treemerger_rf" in node.features:
        faces.add_face_to_node(faces.AttrFace("treemerger_rf", fsize=8), node, 0, position="branch-bottom")

    support_radius= (1.0 - node.support) * 30
    if support_radius > 1:
        support_face = faces.CircleFace(support_radius, "red")
        faces.add_face_to_node(support_face, node, 0, position="float-behind")
        support_face.opacity = 0.25
        faces.add_face_to_node(faces.AttrFace("support", fsize=8), node, 0, position="branch-bottom")
        

    if "clean_alg_mean_identn" in node.features:
        identity = node.clean_alg_mean_identn
    elif "alg_mean_identn" in node.features:
        identity = node.alg_mean_identn

    if "highlighted" in node.features:
        node.img_style["bgcolor"] = "LightCyan"

    if "improve" in node.features:
        color = "orange" if float(node.improve) < 0 else "green"
        if float(node.improve) == 0:
            color = "blue"
             
        support_face = faces.CircleFace(200, color)        
        faces.add_face_to_node(support_face, node, 0, position="float-behind")
        


    
def ncbi_layout(node):
    npr_layout(node)
    global name2color
    if node.is_leaf():
        tax_pos = 10
        if hasattr(node, "lineage"):
            for tax,k in zip(node.lineage, node.named_lineage):
                f = faces.TextFace("%10s" %k, fsize=15)
                try:
                    color = name2color[k]
                except KeyError:
                    name2color[k] = color = treeview.main.random_color()

                #if hasattr(node, "broken_groups") and tax in node.broken_groups:
                f.background.color = color
                faces.add_face_to_node(f, node, tax_pos, position="aligned")
                tax_pos += 1
                    
            f = faces.AttrFace("spname", fsize=15)
            faces.add_face_to_node(f, node, 10, position="branch-right")
    else:
        if getattr(node, "broken_groups", None):
            f = faces.TextFace("\n".join(node.broken_groups), fsize=15, fgcolor="red")
            faces.add_face_to_node(f, node, 1, position="branch-right")

    if hasattr(node, "changed"):
        if node.changed == "yes":
            node.img_style["bgcolor"]="indianred"
        else:
            node.img_style["bgcolor"]="white"

            

def analyze_tracks(t, n2content):
    counterdict = lambda: defaultdict(int)
    node2track = defaultdict(counterdict)
    taxcounter = defaultdict(int)
    tax2name = {}
    for node, leaves in n2content.iteritems():
        if node.is_leaf():
            for index, tax in enumerate(node.lineage):
                taxcounter[tax] += 1
                tax2name[tax] = node.named_lineage[index]
        else:
            for lf in leaves:
                for index, tax in enumerate(lf.lineage):
                    node2track[node][tax] += 1

    mono = set(taxcounter.keys())
    non_mono = set()
    for node, taxa in node2track.iteritems():
        for tax, num in taxa.iteritems():
            if taxcounter[tax] != num and len(n2content[node]) != num:
                if 0:
                    print "max:", taxcounter[tax]
                    print "in this node:", num
                    print "node size:", len(n2content[node])
                    print tax2name[tax]
                    print "..."
                    raw_input()
                mono.discard(tax)
                non_mono.add(tax)
    return mono, non_mono, tax2name
    

def analyze_subtrees(t, subtrees, reft=None):
    ncbi_mistakes = 0
    valid_subtrees = 0
    broken_groups = set()
    correct_groups = set()
    broken_subtrees = 0
    total_rf = 0
    for count, subt in enumerate(subtrees):
        print "\r", count, "   ",
        sys.stdout.flush()
        n2content = subt.get_node2content()
        subt_size = len(n2content[subt])
        if subt_size > 1:
            valid_subtrees += 1
            if reft:
                for _n in n2content():
                    if _n.is_leaf():
                        _n.spcode = _n.realname
                rf, rf_max = subt.robinson_foulds(reft, attr_t1="spcode")
                total_rf += float(rf)/rf_max
            
            si, no, tax2name  = analyze_tracks(subt, n2content)
            ncbi_mistakes += len(no)
            if no:
                broken_subtrees += 1
            broken_groups.update(no)
            correct_groups.update(si)
            children = []
            if args.show_tree or args.render:
                for tip in subt.iter_leaves():
                    target = (t&tip.name)
                    children.append(target)
                    target.broken_groups = set(no)
                # Annotate node
                source_node = t.get_common_ancestor(children)
                source_node.broken_groups = set([tax2name[e] for e in no])
    print "\nDone"
    return valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf

def annotate_tree_with_taxa(t, name2taxa_file, tax2name=None, tax2track=None):
    if name2taxa_file: 
        names2taxid = dict([map(strip, line.split("\t"))
                            for line in open(name2taxa_file)])
    else:
        names2taxid = dict([(n.name, n.name) for n in t.iter_leaves()])

    not_found = 0
    for n in t.iter_leaves():
        n.add_features(taxid=names2taxid.get(n.name, 1))
        n.add_features(species=n.taxid)
        if n.taxid == 1:
            not_found += 1
    if not_found:
        print "WARNING: %s nodes where not found within NCBI taxonomy!!" %not_found
       
    return ncbi.annotate_tree(t, tax2name, tax2track)

    
def tree_compare(t1, t2):
    t2_c2node = {}
    for n, content in t2.get_node2content().iteritems():
        t2_c2node[frozenset([_c.name for _c in content])] = n

    for n, content in t1.get_node2content().iteritems():
        named_content = frozenset([_c.name for _c in content])
        if frozenset(named_content) not in t2_c2node:
            n.add_feature("changed", "yes")
        else:
            n.add_feature("changed", "no")

            
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

    parser.add_argument("--dump", dest="dump",
                        action="store_true", 
                        help="""Dump annotated tree""")

    parser.add_argument("--explore", dest="explore",
                        type=str,
                        help="""Reads a previously analyzed tree and visualize it""")
    
    parser.add_argument("-t", "--tree", dest="target_tree",  nargs="+",
                        type=str, 
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

    parser.add_argument("--rf-only", dest="rf_only",
                        action = "store_true",
                        help="Skip ncbi consensus analysis")

    parser.add_argument("--is_sptree", dest="is_sptree",
                        action = "store_true",
                        help="Assumes no duplication nodes in the tree")
    
    parser.add_argument("-o", dest="output", type=str,
                        help="Writes result into a file")

    parser.add_argument("--tax2name", dest="tax2name", type=str,
                        help="")
    parser.add_argument("--tax2track", dest="tax2track", type=str,
                        help="")
    parser.add_argument("--dump_tax_info", dest="dump_tax_info", action="store_true",
                        help="")
    
    args = parser.parse_args()

    reftree_name = os.path.basename(args.ref_tree) if args.ref_tree else ""
    if args.explore:
        print "Reading tree from file:", args.explore
        t = cPickle.load(open(args.explore))
        ts = TreeStyle()
        ts.force_topology = True
        ts.show_leaf_name = False
        ts.layout_fn = ncbi_layout 
        ts.mode = "r"
        t.show(tree_style=ts)
        print "dumping color config"
        cPickle.dump(name2color, open("ncbi_colors.pkl", "w"))
        sys.exit()
    
    if args.output:
        OUT = open(args.output, "w")
    else:
        OUT = sys.stdout

    print "Dumping results into", OUT
    header = ("Tree".center(50), "Total subtrees", "Broken subtrees", "Broken NCBI clades", "RF (avg)", "RF (med)", "RF (std)", "RF (max possible)")
    print >>OUT, "#"+' '.join([h.center(15) for h in header])
    target_trees = []
    if args.tree_list_file:
        target_trees = [line.strip() for line in open(args.tree_list_file)]
    if args.target_tree:
        target_trees += args.target_tree
    prev_tree = None
    if args.tax2name:
        tax2name = cPickle.load(open(args.tax2name))
    else:
        tax2name = {}

    if args.tax2track:
        tax2track = cPickle.load(open(args.tax2track))
    else:
        tax2track = {}
    print len(tax2track), len(tax2name)
    for tfile in target_trees:
        print tfile
        t = PhyloTree(tfile, sp_naming_function=None)
        t.sort_descendants()

        if prev_tree:
            tree_compare(t, prev_tree)
        prev_tree = t
        
        if args.ref_tree:
            print "Reading ref tree from", args.ref_tree
            reft = Tree(args.ref_tree, format=1)

        else:
            reft = None
        
        if args.tax_info:
            tax2name, tax2track = annotate_tree_with_taxa(t, args.tax_info, tax2name, tax2track)
            if args.dump_tax_info:
                cPickle.dump(tax2track, open("tax2track.pkl", "w"))
                cPickle.dump(tax2name, open("tax2name.pkl", "w"))
                print "Tax info written into pickle files"
        else:
            for n in t.iter_leaves():
                spcode = n.name
                n.add_features(taxid=spcode)
                n.add_features(species=spcode)
            tax2name, tax2track = annotate_tree_with_taxa(t, None, tax2name, tax2track)
                
        # Split tree into species trees
        #subtrees =  t.get_speciation_trees()
        if not args.rf_only:
            print "Calculating tree subparts..."
            t1 = time.time()
            if not args.is_sptree:
                subtrees =  t.split_by_dups()
                print "Subparts:", len(subtrees), time.time()-t1
            else:
                subtrees = [t]
            valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf = analyze_subtrees(t, subtrees)
            print             valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf
        else:
            subtrees = []
            valid_subtrees, broken_subtrees, ncbi_mistakes, total_rf = 0, 0, 0, 0
            
        ndups = 0
        nsubtrees = len(subtrees)
           
        rf = 0
        rf_max = 0
        rf_std = 0
        rf_med = 0
        common_names = 0
        max_size = 0
        if reft and len(subtrees) == 1:
            rf = t.robinson_foulds(reft, attr_t1="realname")
            rf_max = rf[1]
            rf = rf[0]
            rf_med = rf
            
        elif reft:
            print "Calculating avg RF..."
            nsubtrees, ndups, subtrees = t.get_speciation_trees(map_features=["taxid"])
            #print len(subtrees), "Sub-Species-trees found"
            avg_rf = []
            rf_max = 0.0 # reft.robinson_foulds(reft)[1]
            sum_size = 0.0
            print nsubtrees, "subtrees", ndups, "duplications"

            for ii, subt in enumerate(subtrees):
                print "\r%d" %ii,
                sys.stdout.flush()
                try:
                    partial_rf = subt.robinson_foulds(reft, attr_t1="taxid")
                except ValueError:
                    pass
                else:
                    sptree_size = len(set([n.taxid for n in subt.iter_leaves()]))
                    sum_size += sptree_size
                    avg_rf.append((partial_rf[0]/float(partial_rf[1])) * sptree_size)
                    common_names = max(len(partial_rf[2] & partial_rf[3]), common_names)
                    max_size = max(max_size, sptree_size)
                    rf_max = max(rf_max, partial_rf[1])
                #print  partial_rf[:2]
            rf = numpy.sum(avg_rf) / float(sum_size) # Treeko dist
            rf_std = numpy.std(avg_rf)
            rf_med = numpy.median(avg_rf)

        iter_values = (os.path.basename(tfile), reftree_name, nsubtrees, ndups, broken_subtrees, ncbi_mistakes, rf, rf_med, rf_std, rf_max, common_names)
        print >>OUT, '\t'.join(map(str, iter_values))
        OUT.flush()
        if args.show_tree or args.render:
            ts = TreeStyle()
            ts.force_topology = True
            #ts.tree_width = 500
            ts.show_leaf_name = False
            ts.layout_fn = ncbi_layout 
            ts.mode = "r"
            if args.show_tree:
                print "Showing tree..."
                t.show(tree_style=ts)
            else:
                t.render("img.svg", tree_style=ts, dpi=300)
            print "dumping color config"
            cPickle.dump(name2color, open("ncbi_colors.pkl", "w"))

        if args.dump:
            cPickle.dump(t, open("ncbi_analysis.pkl", "w"))
                
    if args.output:
        OUT.close()

    

        
