import argparse
from collections import defaultdict
import sys
import heapq

#############################################################################
# Purpose: find shortest cycle back to (-) node of each (+) node. Start point is (+) node for each (a) strand
# Input: all nodes parsed and renamed from fastg, edges, and weights
# Output: table of all shortest cycles, a separate adjacency matrix for each, and a weight list for r graphs
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_edges', metavar='<list of all renamed edges from fastg>', type=str, help='input fastq files')
parser.add_argument('ifn_edge_summary', metavar='<summary stats for each edge>', type=str, help='input fastq files')
parser.add_argument('out_dir', metavar='<table of mapped k-reads to kmers>', type=str, help='input fastq files')
parser.add_argument('ofn_table', metavar='<adjacency matrix for edge>', type=str, help='input fastq files')

args = parser.parse_args()


def get_out_cov(path,parentToChildren,childToParents):
    sumOutCov = sumInCov = 0
    for node in path:
        outNodes = parentToChildren[node]
        inNodes = childToParents[node]
        for outNode in outNodes:
            if outNode in path:
                continue
            else:
                sumOutCov += edgeToCov[(node,outNode)]
        for inNode in inNodes:
            if inNode in path:
                continue
            else:
                sumInCov += edgeToCov[(inNode,node)]
    return (sumOutCov,sumInCov)

def get_sorted_comp_path(path):
    comp_orientation = {'a+':'b-','b-':'a+','a-':'b+','b+':'a-'}
    comp_path = []
    for node in path:
        contig = node[:-2]
        comp = comp_orientation[node[-2:]]
        comp_path.append(contig+comp)
    return sorted(comp_path)
        
    
print("Collecting all edges now...")

#get all the edges in a hashtable
#find all the starting reference nodes
parentToChildren = {}
childToParents = {}
startNodes = []
with open(args.ifn_edges) as nodes:
    for edge_list in nodes:
        edges = edge_list.rstrip().split(':')
        parent = edges[0]
        children = edges[1].split(',')
        for child in children:
            if child == '':
                continue
            parentToChildren.setdefault(parent,[]).append(child)
            childToParents.setdefault(child,[]).append(parent)
        if parent[-2:] == 'a+':
            startNodes.append(parent)

print("Collecting all edge statistics now...")

edgeToWeight = {}
edgeToCov = {}
edgeToLength = {}
with open(args.ifn_edge_summary) as edge_summary:
    next(edge_summary)
    for edge in edge_summary:
        edgeStats = edge.rstrip().replace(" ","\t").split("\t")
        parent,child,coverage,length,weight = edgeStats
        edge = (parent,child)
        edgeToWeight[edge] = float(weight)
        edgeToCov[edge] = float(coverage)
        edgeToLength[edge] = int(length)

osummary = open(args.ofn_table,'w+')
osummary.write('cycle\ttotal_out_cov\ttotal_in_cov\tavg_cov\tbottleneck_cov\tcycle_length\tedge\tedge_count\tedge_weight\tedge_cov\tedge_length\tnum_contigs\n')
cycleCount = 0
sortedCycleIDs = {}

print("Number of startNodes:",len(startNodes))
startNode_count = 0
for startNode in startNodes:
    startNode_count += 1
    distFromStart = { startNode:0 } #hashtable for distance to start
    unexplored = [ (0,startNode) ] #priority-queue-like heap for unexplored nodes
    explored = {} #keep track of what's been explored
    parentsShortestPath = {} #parent pointers for shortest path.
                             #if the path is relaxed to child, update this pointer
    targetNode = startNode[:-1] + '-' # target node is the parent of the node's invariable edge
    minCycleExists = False
    iteration_count = 0
    queued = {}
    while len(unexplored) > 0: #keep going until there's nothing left to explore
        iteration_count += 1
        minNodeTuple = heapq.heappop(unexplored) #extract min-distance node
        parentNodeDistance,minNode = minNodeTuple
        #print("Node:",startNode,startNode_count,"Iteration:",iteration_count,"Min-node:",minNode,"Dist:",parentNodeDistance,"Visited:",explored.get(minNode))
        if minNode == targetNode: #shortest path to target node reached!
            minCycleExists = True
            break
        explored[minNode] = 'yes' 
        minNodeChildren = parentToChildren.get(minNode)
        if minNodeChildren != None: # if there are children to the latest min-node 
            for child in minNodeChildren: # all children of node being explored
                updatedDistance = False
                childNodeDistance = distFromStart.get(child,sys.maxsize)
                weight = edgeToWeight[(minNode,child)]
                if childNodeDistance > ( parentNodeDistance + weight ): # relax if necessary
                    distance = parentNodeDistance + weight
                    distFromStart[child] = distance
                    parentsShortestPath[child] = minNode # update parent pointer to min distance parent
                    updatedDistance = True
                if explored.get(child) != 'yes': # this child has not been explored (added to the queue)
                    if queued.get(child) == None: # this child has a) not been explored and b) not ever been in the queue
                        queued[child] = 'queued'
                        heapq.heappush( unexplored,( distFromStart[child] , child) ) # only enqueue if hasn't been explored yet and not in queue
                    elif updatedDistance: # this child is already somewhere in the queue, so we have to update its distance (and not add a new instance of it)
                        childIndex = [ind for ind,childTuple in enumerate(unexplored) if childTuple[1] == child] # find the index of the child tuple so we cant delete it
                        unexplored[childIndex[0]] = (distFromStart[child],child)  # update heap element of child to potential new distance from start
                        heapq.heapify(unexplored) # maintain heap after updating child distance
                        
    if minCycleExists: #if there's a min weight cycle
        shortestPath = [targetNode]
        childNode = targetNode
        while childNode != startNode:
            parentNode = parentsShortestPath[childNode]
            shortestPath.append(parentNode)
            childNode = parentNode
        sortedID = tuple(sorted(shortestPath))
        if sortedCycleIDs.get(sortedID) == 'added': # check if the cycle already exists, skip if it does
            continue
        sortedCycleIDs[sortedID] = 'added' # make sure that future duplicate cycles aren't returned
        sortedCycleIDs[tuple(get_sorted_comp_path(sortedID))] = 'added' # also add the complement path, which is reverse complement
        cycleCount += 1
        
        shortestPath.append(targetNode)
        shortestPath.reverse()
    
        (sumOutCov,sumInCov) = get_out_cov(shortestPath[:-1],parentToChildren,childToParents)
        cycleLength = 0
        cycleBottleNeck = sys.maxsize
        numEdges = 0
        cycleCoverage = 0
        for parentIndex in range(0,len(shortestPath)-1):
            numEdges += 1
            edge = (shortestPath[parentIndex],shortestPath[parentIndex+1])
            edgeLength = edgeToLength[edge]
            cycleLength += edgeLength
            edgeCoverage = edgeToCov[edge]
            cycleCoverage += edgeCoverage
            if edgeCoverage < cycleBottleNeck:
                cycleBottleNeck = edgeCoverage
        avgCycleCov = float(cycleCoverage)/numEdges
        edge_count = 0
        cycleLabel = 'CYC' + str(cycleCount)
        numContigsInCyc = str((len(shortestPath)-1) / 2)
        for parentIndex in range(0,len(shortestPath)-1):
            edge_count += 1
            edge = (shortestPath[parentIndex],shortestPath[parentIndex+1])
            edgeLabel = edge[0] + '>' + edge[1]
            edgeLength = edgeToLength[edge]
            edgeCoverage = edgeToCov[edge]
            edgeWeight = edgeToWeight[edge]
            
            osummary.write(cycleLabel + "\t" + str(sumOutCov) + "\t" + str(sumInCov) + "\t" + str(avgCycleCov) + "\t" + str(cycleBottleNeck) + '\t' +
                            str(cycleLength) + "\t" + edgeLabel + "\t" + str(edge_count) + "\t" +
                            str(edgeWeight) + "\t" + str(edgeCoverage) + "\t" + str(edgeLength) + "\t" + numContigsInCyc + "\n")
        
        omatrix = open(args.out_dir + "/" + cycleLabel + "_adjmatrix.txt", "w+")
        oedgeweights = open(args.out_dir + "/" + cycleLabel + "_edgeweights.txt", "w+")
        omatrix.write("\t")
        for node in shortestPath[:-1]:
            omatrix.write(node + "\t")
        omatrix.write("\n")
        for parent in shortestPath[:-1]:
            omatrix.write(parent + "\t")
            children = parentToChildren[parent]
            for child in shortestPath[:-1]:
                if child in children:
                    connection = "1"
                    edge = (parent,child)
                    weight = edgeToWeight[edge]
                    cov = edgeToCov[edge]
                    oedgeweights.write(edge[0] + ">" + edge[1] + "\t" + "W:" + str(int(weight)) + ",C:" + str(int(cov)) + "\n")
                else:   
                    connection = "0"
                omatrix.write(connection + "\t")
            omatrix.write("\n")
            
            
            
    
    
            
            
            
            
        
    
        
