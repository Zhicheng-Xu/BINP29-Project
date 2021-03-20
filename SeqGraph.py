#!/usr/bin/python3
"""
 Title: SeqGraph.py
 Date: 2021-03-16
 Authors: Zhicheng Xu
 
 Description:
     This program will split the input sequence into K-mers for input K value
     and generate a De Bruijn graph and reassemble the sequence to calculate
     the similarity between the original sequence and the new one.
     
 Usage:
     ./alignment.py
"""
import re
from flask import Flask, render_template, request, url_for
from graphviz import Digraph
    
app = Flask(__name__)


class Node:
    def __init__(self, label):
        self.label = label
        self.indeg = 0
        self.outdeg = 0

def isEulerian(nodes): #check if the graph has a Eulerian path
    semibalanced = 0
    balanced = 0
    for label in nodes.keys():
        if nodes[label].indeg == nodes[label].outdeg:
            balanced += 1
        elif nodes[label].outdeg - nodes[label].indeg == 1:
            semibalanced += 1
            start = nodes[label]
        elif nodes[label].indeg - nodes[label].outdeg == 1:
            semibalanced += 1
            end = nodes[label]
    
    if semibalanced == 2 and semibalanced+balanced == len(nodes):
        return(True, start, end)
    else:
        return(False, None, None)
 


@app.route('/')
def sequence_input():
    return render_template('sequence_input.html', 
                           bg_src = url_for('static', filename='network.png'),
                           msg = "")

@app.route('/', methods=['POST', 'GET'])
def reassemble():
    path = []
    def findPath(node):
        while new_graph[node]: #stop until the list is empty
            toNode = new_graph[node].pop(0) #walk through the list and remove the visited node
            findPath(toNode) #recursive method
        path.append(node)
    
    if request.method == 'POST':
        seq = request.form['seq']
        seq = seq.strip()
        if not seq:
            message = 'Please input valid sequence and K value!'
            return render_template('sequence_input.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   msg = message)
        
        value = re.search('[^-ATCGatcg]', seq) 
        if value != None: #check if the input sequence is valid
            message = 'Your sequence is not the pure nucleotide!'
            return render_template('sequence_input.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   msg = message)
    
        k = request.form['Kmer']
        k = k.strip()
        if not k: #check if the input K value is valid
            message = 'Please input valid sequence and K value!'
            return render_template('sequence_input.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   msg = message)
        
        if not k.isnumeric():
            message = 'Your K value is not an integer!'
            return render_template('sequence_input.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   msg = message)
        
        k = int(k)
        if k > len(seq):
            message = 'Your K value is too large!'
            return render_template('sequence_input.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   msg = message)
        
        graph = {}
        nodes = {}
        
        seq = seq.upper()
        for i in range(len(seq)-k+1): #for each K-mer
            label1 = seq[i:i+k-1] #left K-1 mer
            label2 = seq[i+1:i+k] #right K-1 mer
            
            if label1 not in nodes.keys(): #put the nodes into a set
                node1 = Node(label1)       
                nodes[label1] = node1
            else:
                node1 = nodes[label1]
                
            if label2 not in nodes.keys():
                node2 = Node(label2)
                nodes[label2] = node2
            else:
                node2 = nodes[label2]
                
            node1.outdeg += 1
            node2.indeg += 1
            
            if node1 not in graph.keys(): #put the edge into a dictionary
                graph[node1] = [node2]  #start node as key, end node as value of a list
            else:
                graph[node1].append(node2)
            
                    
        dot = Digraph(format = 'png') #visualize the graph
        for vertex in nodes.keys():
            dot.node(vertex)
            
        for key in graph.keys():
            for item in graph[key]:
                dot.edge(key.label, item.label)
        
        img_src = 'digraph_'+str(k)+'.gv'
        dot.render('static/'+img_src, view=False)
        
        
        if not isEulerian(nodes)[0]: #check if the graph is Eulerian
            message = "The De Bruijn graph is not eulerian. Please choose another K value for Kmers!"
            return render_template('sequence_input.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   msg = message)
        
        else:
            start = isEulerian(nodes)[1]
            end = isEulerian(nodes)[2]
            
            new_graph = graph.copy()
            if end not in new_graph.keys(): #add an edge from the end Node to the start Node
                new_graph[end] = [start]
            else:
                new_graph[end].append(start)
                
    
            findPath(start) #find the path recursively
            path = path[::-1] #reverse the path and remove the added edge
            new_seq = start.label
            for node in path[1:-1]:
                new_seq += node.label[-1]
                
            similarity = 0
            if len(seq) == len(new_seq): #calculate the similarity between the old sequence and the new one
                for i in range(len(seq)):
                    if seq[i] == new_seq[i]:
                        similarity += 1
                similarity = similarity/len(seq)*100
            else:
                if len(seq) > len(new_seq):
                    for i in range(len(new_seq)):
                        if seq[i] == new_seq[i]:
                            similarity += 1
                    similarity = similarity/len(seq)*100
                else:
                    for i in range(len(seq)):
                        if seq[i] == new_seq[i]:
                            similarity += 1
                    similarity = similarity/len(new_seq)*100
            
            return render_template('show_graph.html', 
                                   bg_src = url_for('static', filename='network.png'),
                                   img_src = url_for('static', filename=img_src+'.png'), 
                                   old_seq = seq,
                                   new_seq = new_seq,
                                   similarity = similarity) 
    else:
        message = 'Please input valid sequence and K value!'
        return render_template('sequence_input.html', msg = message)
    
if __name__ == '__main__':
    app.run(debug=True)
    
    
    
    
    
    
    
    
    