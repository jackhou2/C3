# -*- coding: utf-8 -*-

"""

Created on Tue Nov 11 10:26:46 2014



@author: Jack Hou and Amin Emad

"""



from __future__ import division



import community_detection_EMP_PAN as cd

import numpy as np

import networkx as nx

#import scipy as sp

#import matplotlib as mpl

#import matplotlib.pyplot as plt

import itertools as its

import pandas as pd

import csv

import copy

import time

import pickle



############################################################

def ts(*x):      #the *x says that several variables may be passed to the function

    return tuple(sorted(x))



############################################################

def CSV_imp(address):

    with open(address,"r") as fin:

        data_reader = list(csv.reader(fin,delimiter=','))

        data_reader.pop(0)

    return (data_reader)     









############################################################

def Gene_filter(address_mutation,address_copynumber,t_percentile=90,t_cnv_up=0,t_cnv_low=0):

    mutation_array = (np.genfromtxt(address_mutation,delimiter=','))[1:,1:]

    copynumber_array = (np.genfromtxt(address_copynumber,delimiter=','))[1:,1:]

    y = np.logical_or((mutation_array>0.5),np.logical_or((copynumber_array>t_cnv_up+0.5),(copynumber_array<t_cnv_low-0.5)))

    temp = y.sum(axis=1)

    t_filter = np.percentile(temp,t_percentile)

    y_temp = np.arange(0,len(temp))

    y_label = y_temp[temp >= t_filter]

    y_pruned = np.array(y[y_label][:],'int8')

    #Extracting gene names

    with open(address_mutation,"r") as fin:

        data_reader = list(csv.reader(fin,delimiter=','))

        data_reader.pop(0)        

    gene_names = [row[0] for row in data_reader]    

    

    return (gene_names,y_label,y_pruned,t_filter)



############################################################

def gene_state(y):

    gpair_state = {}   #

    """This is a dictionary of the form (i,j): (n11,n10,min(ni,nj))"""        

    numgene = np.size(y,axis=0)

    maxval = 0

    for (i,j) in its.combinations(range(numgene),2):  

            gpair_state[ts(i,j)] = (((y[i,]*y[j,]).sum(axis=0)),((y[j,]+y[i,]-2*y[i,]*y[j,]).sum(axis=0)),np.minimum(np.sum(y[j,]),np.sum(y[i,])))

            if gpair_state[ts(i,j)][0]>maxval:

                maxval = copy.copy(gpair_state[ts(i,j)][0])

    return (maxval,gpair_state)     

       

       

############################################################

def weight_CO_ME(y_label,gpair_state,numgene,a,J_percentile):

    #y_label shows the label of remaining genes

    #G2 is the complete graph of remaining genes

    G2 = nx.complete_graph(numgene)

    

    wn = {ts(u,v): gpair_state[ts(u,v)][0] / gpair_state[ts(u,v)][2] for (u,v) in G2.edges()}   

    #wp = {ts(u,v): gpair_state[ts(u,v)][1] for (u,v) in G2.edges()}

    

    

    #wp_max_temp = max([wp[u] for u in wp]) 

    #wn_max_temp = max([wn[u] for u in wn])    

    

    wp = {}

    gpair10_array = np.array([gpair_state[u][1] for u in gpair_state])

    threshold = np.percentile(gpair10_array,J_percentile)

    for x in gpair_state:

        wn[x] = a * wn[x]

        if gpair_state[x][1]>= threshold:

            wp[x] = 1

        else:

            wp[x] = 1/threshold * gpair_state[x][1]

            b = wp[x] + wn[x]

            if b < 1:

                if b > 0:

                    wn[x] = wn[x] / b

                    wp[x] = 1 - wn[x]

                else:

                    wn[x] = 0.5

                    wp[x] = 0.5

                

    nx.set_edge_attributes(G2,"w_n",wn)

    nx.set_edge_attributes(G2,"w_p",wp)

    return (G2,wp,wn) 

    

############################################################

def weight_NI_ME(address_network,y_label,gpair_state,numgene,a,J_percentile):

    #y_label shows the label of pruned genes

    #address_network is the address of the file containing adjecency matrix of the HRN network

    #We assume that the file is in CSV format and has a header and a column sider 

    #G2 is the complete graph of important genes

    G2 = nx.complete_graph(numgene)



    wn = {ts(u,v): gpair_state[ts(u,v)][0] / gpair_state[ts(u,v)][2] for (u,v) in G2.edges()}   

     

    net_adjacency = (np.genfromtxt(address_network,delimiter=','))[1:,1:]  

    wp = {}

    for u in range(numgene):

        y_tempu = np.logical_or(net_adjacency[y_label[u],],net_adjacency[:,y_label[u]])

        #we find the logical or of the y_label[u]'th column and y_label[u]'th row of the adjacency matirx         

        y_tempu[y_label[u]] = True

        for v in range(u+1,numgene):

            y_tempv = np.logical_or(net_adjacency[y_label[v],],net_adjacency[:,y_label[v]])

            y_tempv[y_label[v]] = True

            wp[ts(u,v)] = sum(np.logical_and(y_tempu,y_tempv))/sum(np.logical_or(y_tempu,y_tempv))



    wp_array = np.array([wp[u] for u in wp]) 

    threshold = np.percentile(wp_array,J_percentile)

    #wp_max_temp = max([wp[u] for u in wp]) 

    #wn_max_temp = max([wn[u] for u in wn])    

    

    for x in wp:

        wn[x] = a * wn[x]



        if wp[x]>= threshold:

            wp[x] = 1

        else:

            wp[x] = 1/threshold * wp[x]

            b = wp[x] + wn[x]

            if b < 1:

                if b > 0:

                    wn[x] = wn[x] / b

                    wp[x] = 1 - wn[x]

                else:

                    wn[x] = 0.5

                    wp[x] = 0.5

    nx.set_edge_attributes(G2,"w_n",wn)

    nx.set_edge_attributes(G2,"w_p",wp)

    return (G2,wp,wn) 

############################################################

def weight_NI_CO_ME(address_network,y_label,gpair_state,numgene,a,J_percentile,c):

    #y_label shows the label of pruned genes

    #address_network is the address of the file containing adjecency matrix of the HRN network

    #We assume that the file is in CSV format and has a header and a column sider 

    #G2 is the complete graph of important genes

    G2 = nx.complete_graph(numgene)

    

    wn = {ts(u,v): gpair_state[ts(u,v)][0] / gpair_state[ts(u,v)][2] for (u,v) in G2.edges()}       

    

    wp1 = {}

    gpair10_array = np.array([gpair_state[u][1] for u in gpair_state])

    threshold = np.percentile(gpair10_array,J_percentile)

    for x in gpair_state:

        if gpair_state[x][1]>= threshold:

            wp1[x] = 1

        else:

            wp1[x] = 1/threshold * gpair_state[x][1] 

        

    net_adjacency = (np.genfromtxt(address_network,delimiter=','))[1:,1:]    

    wp2 = {}

    for u in range(numgene):

        y_tempu = np.logical_or(net_adjacency[y_label[u],],net_adjacency[:,y_label[u]])

        #we find the logical or of the y_label[u]'th column and y_label[u]'th row of the adjacency matirx         

        y_tempu[y_label[u]] = True

        for v in range(u+1,numgene):

            y_tempv = np.logical_or(net_adjacency[y_label[v],],net_adjacency[:,y_label[v]])

            y_tempv[y_label[v]] = True

            wp2[ts(u,v)] = sum(np.logical_and(y_tempu,y_tempv))/sum(np.logical_or(y_tempu,y_tempv))

            

    wp2_array = np.array([wp2[u] for u in wp2]) 

    threshold2 = np.percentile(wp2_array,J_percentile)

    

    for x in wp2:

        if wp2[x]>= threshold2:

            wp2[x] = 1

        else:

            wp2[x] = 1/threshold2 * wp2[x]

    wp = {ts(u,v): c*(wp1[ts(u,v)]) + (1-c)*(wp2[ts(u,v)]) for (u,v) in G2.edges()}

    for x in wp:

        wn[x] = a * wn[x]

        b = wp[x] + wn[x]

        if b < 1:

            if b > 0:

                wn[x] = wn[x] / b

                wp[x] = 1 - wn[x]

            else:

                wn[x] = 0.5

                wp[x] = 0.5

       

    nx.set_edge_attributes(G2,"w_n",wn)

    nx.set_edge_attributes(G2,"w_p",wp)

    return (G2,wp,wn) 

    

############################################################

def weight_EX_CO_ME(address_expression,y_label,gpair_state,numgene,a,J_percentile,c):

    #y_label shows the label of pruned genes

    #address_expression is address of the file containing z-score of the expressions

    #We assume that the file is in CSV format and has a header and a column sider 

    #G2 is the complete graph of important genes
    G2 = nx.complete_graph(numgene)

    wn = {ts(u,v): gpair_state[ts(u,v)][0] / gpair_state[ts(u,v)][2] for (u,v) in G2.edges()}       
    wp1 = {}
    gpair10_array = np.array([gpair_state[u][1] for u in gpair_state])

    threshold = np.percentile(gpair10_array,J_percentile)

    for x in gpair_state:

        if gpair_state[x][1]>= threshold:

            wp1[x] = 1

        else:

            wp1[x] = 1/threshold * gpair_state[x][1] 

        

    expression = (np.genfromtxt(address_expression,delimiter=','))[1:,1:] 

    

    wp2 = {}

    for (u,v) in G2.edges():

        temp_cos = np.absolute(np.dot(expression[y_label[u],],expression[y_label[v],])) / np.sqrt(np.dot(expression[y_label[u],],expression[y_label[u],]) * np.dot(expression[y_label[v],],expression[y_label[v],])) 

        if np.isnan(temp_cos):

            wp2[ts(u,v)] = 0.5

        else:

            wp2[ts(u,v)] = temp_cos

            

    wp2_array = np.array([wp2[u] for u in wp2]) 

    threshold2 = np.percentile(wp2_array,J_percentile)

    

    for x in wp2:

        if wp2[x]>= threshold2:

            wp2[x] = 1

        else:

            wp2[x] = 1/threshold2 * wp2[x]

             

    wp = {ts(u,v): c*(wp1[ts(u,v)]) + (1-c)*(wp2[ts(u,v)]) for (u,v) in G2.edges()}

    

    for x in wp:

        wn[x] = a * wn[x]



        b = wp[x] + wn[x]

        if b < 1:

            if b > 0:

                wn[x] = wn[x] / b

                wp[x] = 1 - wn[x]

            else:

                wn[x] = 0.5

                wp[x] = 0.5

       

    nx.set_edge_attributes(G2,"w_n",wn)

    nx.set_edge_attributes(G2,"w_p",wp)

    return (G2,wp,wn) 

    

############################################################

def gene_map(gene_names,rounded,y_label):

    final_cluster = []  

    for row in rounded:

        temp = [gene_names[y_label[u]] for u in row]

        final_cluster.append(temp)

    return final_cluster


#############################################################




def weight_Triple(address_expression,address_network,y_label,gpair_state,numgene,a,J_percentile,cov,exp,net):
    G2 = nx.complete_graph(numgene)
    wn = {ts(u,v): gpair_state[ts(u,v)][0] / gpair_state[ts(u,v)][2] for (u,v) in G2.edges()}       

    #convert all weights to a 0-1 number
    sumweight = cov + exp + net
    cov = cov/sumweight
    exp = exp/sumweight
    net = net/sumweight
    print cov
    print exp
    print net


#This is the Wp1 state

    wp1 = {}
    gpair10_array = np.array([gpair_state[u][1] for u in gpair_state])
    threshold = np.percentile(gpair10_array,J_percentile)
    for x in gpair_state:
        if gpair_state[x][1]>= threshold:
            wp1[x] = 1
        else:
            wp1[x] = 1/threshold * gpair_state[x][1] 


#This is the Wp2 state

    expression = (np.genfromtxt(address_expression,delimiter=','))[1:,1:] 
    wp2 = {}

    for (u,v) in G2.edges():
        temp_cos = np.absolute(np.dot(expression[y_label[u],],expression[y_label[v],])) / np.sqrt(np.dot(expression[y_label[u],],expression[y_label[u],]) * np.dot(expression[y_label[v],],expression[y_label[v],])) 
        if np.isnan(temp_cos):
            wp2[ts(u,v)] = 0.33
        else:
            wp2[ts(u,v)] = temp_cos

    wp2_array = np.array([wp2[u] for u in wp2]) 
    threshold2 = np.percentile(wp2_array,J_percentile)

    for x in wp2:
        if wp2[x]>= threshold2:
            wp2[x] = 1
        else:
            wp2[x] = 1/threshold2 * wp2[x]


#This is the Wp3 state
    net_adjacency = (np.genfromtxt(address_network,delimiter=','))[1:,1:]  

    wp3 = {}

    for u in range(numgene):

        y_tempu = np.logical_or(net_adjacency[y_label[u],],net_adjacency[:,y_label[u]])

        #we find the logical or of the y_label[u]'th column and y_label[u]'th row of the adjacency matirx         

        y_tempu[y_label[u]] = True

        for v in range(u+1,numgene):

            y_tempv = np.logical_or(net_adjacency[y_label[v],],net_adjacency[:,y_label[v]])

            y_tempv[y_label[v]] = True

            wp3[ts(u,v)] = sum(np.logical_and(y_tempu,y_tempv))/sum(np.logical_or(y_tempu,y_tempv))

    wp3_array = np.array([wp3[u] for u in wp3]) 

    threshold3 = np.percentile(wp3_array,J_percentile)

    for x in wp3:

        if wp3[x]>= threshold3:

            wp3[x] = 1

        else:

            wp3[x] = 1/threshold3 * wp3[x]





    wp = {ts(u,v): cov*(wp1[ts(u,v)]) + exp*(wp2[ts(u,v)]) + net*(wp3[ts(u,v)]) for (u,v) in G2.edges()}
    #w1 = cov*(wp1[ts(u,v)])
    #w2 = exp*(wp2[ts(u,v)])
    #w3 = net*(wp3[ts(u,v)])

    #print w1
    #print w2
    #print w3

    for x in wp:
        wn[x] = a * wn[x]
        b = wp[x] + wn[x]
        if b < 1:
            if b > 0:
                wn[x] = wn[x] / b
                wp[x] = 1 - wn[x]
            else:
                wn[x] = 0.5
                wp[x] = 0.5
    
    nx.set_edge_attributes(G2,"w_n",wn)
    nx.set_edge_attributes(G2,"w_p",wp)
    return (G2,wp,wn)




###################MAIN

start_time = time.time()



###################MAIN
import argparse

parser = argparse.ArgumentParser()

#Add the cancer argument

parser.add_argument("--cancer", help="Add a cancer",
                    action="store")

parser.add_argument("--filter", type=float,help="filter. Default value 20",
                    action="store")

parser.add_argument("--percentile",type=float, help="percentile. Default value 90",
                    action="store")

parser.add_argument("--up", type=int,help="up. Default value 3",
                    action="store")

parser.add_argument("--down",type=int, help="down. Default value -1",
                    action="store")

parser.add_argument("--J_percentile", type=float, help="J_percentile. Default value 90",
                    action="store")

parser.add_argument("--alpha", type=float, help="alpha. Default value 0.29",
                    action="store")

parser.add_argument("--date", help="date. Default value April29",
                    action="store")

parser.add_argument("--a", type=int, help="a. Default value 4",
                    action="store")

parser.add_argument("--bound", type=int, help="bound. Default value 5",
                    action="store")

parser.add_argument("--cond", help="cond. Default value EX_CO_ME. Triple for all three of EX, CO, ME",
                    action="store")

parser.add_argument("--address", help="address. Default value EX_CO_ME",
                    action="store")

parser.add_argument("--output", help="output",
                    action="store")

parser.add_argument("--c",type=float,help="c, the weight of the variable",
                    action="store")

parser.add_argument("--network",help="the network",action="store")


parser.add_argument("--cov", type=float, help="coverage weight. For weight_Triple only. Default 0.33",
                    action="store")

parser.add_argument("--exp", type=float, help="expression weight. For weight_Triple only.For  Default 0.33",
                    action="store")

parser.add_argument("--net", type=float, help="coverage weight. Default 0.33",
                    action="store")




args = parser.parse_args()
if args.cancer:
   cancer = args.cancer
else:
   raise ValueError('no cancer added')

print "cancer is %s" % cancer

#Add the filter argument

if args.filter:
   t_filter = args.filter
else:
   t_filter = 20

print "t_filter is %s" % t_filter


#Add the percentile argument

if args.percentile:
   t_percentile = args.percentile
else:
   t_percentile = 90

print "percentile is %s" % t_percentile

#Add the cnv up argument

if args.up:
   t_cnv_up = args.up
else:
   t_cnv_up = 3

print "t_cnv_up is %s" % t_cnv_up


#Add the cnv down argument

if args.up:
   t_cnv_low = args.down
else:
   t_cnv_low = -1

print "t_cnv_low is %s" % t_cnv_low

#Add the j percentile argument

if args.J_percentile:
   J_percentile = args.J_percentile
else:
   J_percentile = 90

print "J_percentile is %s" % J_percentile

#Add the alpha argument

if args.alpha:
   alpha = args.alpha
else:
   alpha = 0.29

print "alpha is %s" % alpha


#Add the date argument

if args.date:
   date = args.date
else:
   date = "April29"

print "date is %s" % date


#Add the a argument

if args.a:
   a = args.a
else:
   a = 4

print "a is %s" % a


#Add the bound argument

if args.bound:
   bound = args.bound
else:
   bound = 5

print "bound is %s" % bound

#Add the bound argument

if args.cond:
   cond = args.cond
else:
   cond = "EX_CO_ME"

print "cond is %s" % cond


#Add the address argument

if args.address:
   address = args.address
else:
   raise ValueError('no address added')

print "address is %s" % address

#Add the address argument

if args.output:
   output = args.output
else:
   raise ValueError('no output added')

print "output is %s" % output

#Add the weight argument

if args.c:
   c = args.c
else:
   c = 0.5

print "c is %s" % c

if args.network:
   network = args.network
else:
   network = address

print "network is %s" % network


if args.cov:
   cov = args.cov
else:
   cov = 0.33

print "coverage weight is %s" % cov

if args.exp:
   exp = args.exp
else:
   exp = 0.33

print "expression weight is %s" % exp

if args.net:
   net = args.net
else:
   net = 0.33

print "network weight is %s" % net


print "Cancer = %s, t_percentile = %s, t_cnv_up = %s, t_cnv_low = %s, a = %s, size bound = %s, alpha = %s, conditions = %s" %(cancer,t_percentile,t_cnv_up, t_cnv_low,a,bound,alpha,cond)
            
            # Addresses for Prof. Ma server

address_mutation = '%s/%s/%s_mutation.csv' %(address,cancer,cancer)
address_copynumber = '%s/%s/%s_copyNumber.csv' %(address,cancer,cancer)
address_network = "%s/OrderedKEGGNetwork.csv" % network
address_expression = '%s/%s/%s_expression.csv' %(address,cancer,cancer)
pickle_name = "%s/pickle_%s_t%s_tcnvu%s_tcnvl%s_a%s_alpha_%s_bound%s_%s_%s_%s_EMP.pickle" %(output,cancer,t_percentile,t_cnv_up,t_cnv_low,a,alpha,bound,date,cond,c)
            #address_out = "/home/emad2/programs/output_%s_t%s_a%s_bound%s_%s.csv" %(cancer,t_filter,a,bound,date,c)
        

(gene_names,y_label,y_pruned,t_filter) = Gene_filter(address_mutation,address_copynumber,t_percentile,t_cnv_up,t_cnv_low)
            #y = (pd.read_table(address_mutation)).as_matrix()
            #(y_label,y_pruned) = cd.mut_prune(y,t_filter)
            
            	#a = 4

            	#bound = 5

            

print "Cancer = %s, t_percentile = %s, t_cnv_up = %s, t_cnv_low = %s, a = %s, size bound = %s, alpha = %s, conditions = %s" %(cancer,t_percentile,t_cnv_up, t_cnv_low,a,bound,alpha,cond)

            
(gene_names,y_label,y_pruned,t_filter) = Gene_filter(address_mutation,address_copynumber,t_percentile,t_cnv_up,t_cnv_low)

            		#y = (pd.read_table(address_mutation)).as_matrix()

            		#(y_label,y_pruned) = cd.mut_prune(y,t_filter)

numgene = np.size(y_pruned, axis=0)

print "number of genes is", numgene

numpat = np.size(y_pruned, axis = 1)

print "number of samples is", numpat

(maxval,gpair_state) = gene_state(y_pruned)

print time.time()-start_time, "seconds for stage 1"

print address_network

if cond == "NI_ME":

	(G,wp,wn) = weight_NI_ME(address_network,y_label,gpair_state,numgene,a,J_percentile)

elif cond == "CO_ME":

	(G,wp,wn) = weight_CO_ME(y_label,gpair_state,numgene,a,J_percentile)

elif cond == "NI_CO_ME":

	(G,wp,wn) = weight_NI_CO_ME(address_network,y_label,gpair_state,numgene,a,J_percentile,c)

elif cond == "EX_CO_ME":

	(G,wp,wn) = weight_EX_CO_ME(address_expression,y_label,gpair_state,numgene,a,J_percentile,c)

elif cond == "triple":

	(G,wp,wn) = weight_Triple(address_expression,address_network,y_label,gpair_state,numgene,a,J_percentile,cov,exp,net)
	pickle_name = "%s/pickle_%s_t%s_tcnvu%s_tcnvl%s_a%s_alpha_%s_bound%s_%s_%s_COV_%s_EXP_%s_NET_%s_EMP.pickle" %(output,cancer,t_percentile,t_cnv_up,t_cnv_low,a,alpha,bound,date,cond,cov,exp,net)




else: 

	print "ERROR"

print time.time()-start_time, "seconds for stage 2"

(clust, frac_cost) = cd.optimal_clustering(G, size_bound = bound)

print time.time()-start_time, "seconds for stage 3"

(rounded, clust_rounded, rounded_cost) = cd.round_solution(G, numgene, clust,alpha,size_bound = bound)

print time.time()-start_time, "seconds for stage 4"

gene_clusters = gene_map(gene_names,rounded,y_label)

print time.time()-start_time, "seconds for stage 5"

with open(pickle_name,"w") as f:

	pickle.dump([J_percentile,a, G, clust_rounded, rounded_cost, alpha, bound, cancer, clust, frac_cost, gene_clusters, gene_names, gpair_state, maxval, numgene, numpat, rounded, t_filter, t_percentile, t_cnv_up, t_cnv_low, wn, wp, y_label, y_pruned],f)
            

print time.time()-start_time, "seconds total"

            
 

