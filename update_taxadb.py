from string import strip
from ete2a1 import Tree

def load_ncbi_tree_for_dump():
    # Download: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    parent2child = {}
    name2node = {}
    node2taxname = {}
    for line in open("names.dmp"):
        fields =  map(strip, line.split("|"))
        nodename = fields[0]
        name_type = fields[3].lower()
        taxname = fields[1]
        if name_type == "scientific name":
            node2taxname[nodename] = taxname 

    print len(node2taxname), "names loaded"

    for line in open("nodes.dmp"):
        fields =  line.split("|")
        nodename = fields[0].strip()
        parentname = fields[1].strip()
        n = Tree()
        n.name = nodename
        n.taxname = node2taxname[nodename]
        parent2child[nodename] = parentname
        name2node[nodename] = n

    print len(name2node), "nodes loaded"

    for node in name2node:
       if node == "1":
           t = name2node[node]
       else:
           parent = parent2child[node]
           parent_node = name2node[parent]
           parent_node.add_child(name2node[node])
    return t

def generate_table(t):
    OUT = open("taxa.tab", "w")
    for j, n in enumerate(t.traverse()):
        if j%1000 == 0:
            print "\r",j,
        temp_node = n
        track = []
        while temp_node:
            track.append(temp_node.name)
            temp_node = temp_node.up
        if n.up:
            print >>OUT, '\t'.join([n.name, n.up.name, n.taxname, ','.join(track)])
        else:
            print >>OUT, '\t'.join([n.name, "", n.taxname, ','.join(track)])
    OUT.close()

t = load_ncbi_tree_for_dump()
print "printing:",
generate_table(t)
print """ 
 # sqlite3 taxa.sqlite
 create table species
 (taxid INT PRIMARY KEY, parent INT, spname VARCHAR(50), track TEXT);
  .separator "\t"
  .import taxa.tab species
"""