# BINP29-Project
SeqGraph is a software to stimulate the assembly process of Velvet. This program will split the input sequence into K-mers for input K valueand generate a De Bruijn graph and reassemble the sequence to calculate the similarity between the original sequence and the new one.  
Users need to input a nucleotide sequence and an integer for K value to gennerate the corresponding De Bruijn graph and sequence assembly. The results would be displayed on a web page based on web programming.  
The 'static' directory is for image storage. For each different sequence input, users need to clear this directory manually except for the background image 'network.png' which is used to build the web page. This is because browser cache would remember the existing imaged even though the image has been replaced in the directory.  
The 'templates' directory is for templates storage of the web program.  
For any questions, please contact zh2680xu-s@student.lu.se  
  
###### Python Implement
Initialize a class named 'Node' that has three attritubes: label for K-1 mer, indegree for the number of edges that point to this node, and outdegree for the number of edges that point out of this node.  

    class Node:
      def __init__(self, label):
          self.label = label
          self.indeg = 0
          self.outdeg = 0

Define a function to check if the graph is Eulerian.

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
          
Based on web programming, the main page of the program is defined as sequence_input().

    @app.route('/')
    def sequence_input():
        return render_template('sequence_input.html', 
                               bg_src = url_for('static', filename='network.png'),
                               msg = "")

After the user enter the input sequence and K value and click the 'submit' button, the program will first check the validity of the input sequence and K value.

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
                                   
For each sequence, initialize a dictionary called 'graph' that takes the start Node of each edge as the key, and the end Nodes list which the key Node points to as the value. Another dictionary is to map the label to 'Node'.  
For each K-mer, each node represents a distinct K-1 mer and each K-mer has one directed edge from the left K-1 mer to the right K-1 mer.

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
                
To visualize the graph, Graphvia package is useful.

    dot = Digraph(format = 'png') #visualize the graph
        for vertex in nodes.keys():
            dot.node(vertex)
            
        for key in graph.keys():
            for item in graph[key]:
                dot.edge(key.label, item.label)
        
        img_src = 'digraph_'+str(k)+'.gv'
        dot.render('static/'+img_src, view=False)
        
Check if the graph is Eulerian and recursively walk through the graph.

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

    
