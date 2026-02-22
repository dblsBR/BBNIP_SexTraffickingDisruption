import networkx as nx;
import gurobipy as gp;
from gurobipy import GRB;
from datetime import datetime;
import math;
import numpy;
from numpy import * 


def Mod_MFNIP (G, s, t, budget, vic, Gamma2, restriction):     
    
    model = gp.Model("Modified_SpecialArcs"); 
    
    objVar = model.addVar(vtype=GRB.INTEGER, name="objVar", lb = 0, ub = GRB.INFINITY);
    
    gamma = model.addVars(G.edges, vtype=GRB.BINARY, name="gamma"); 
    
    x = model.addVars(G.edges, vtype=GRB.CONTINUOUS, name="x", lb = 0, ub = GRB.INFINITY);
    theta = model.addVars(G.edges, vtype=GRB.CONTINUOUS, name="theta", lb = 0, ub = 1); 
    alpha = model.addVars(G.nodes, vtype=GRB.CONTINUOUS, name="alpha", lb = 0, ub = 1); 
    beta = model.addVars(G.edges, vtype=GRB.CONTINUOUS, name="beta", lb = 0, ub = GRB.INFINITY); 
    
    
    # Restriction on Interdiction
    if restriction == 1:
        for i,j in G.edges:
            if G.edges[i,j]['trafficker'] == 1:
                model.addConstr(gamma[i,j] == 0);
            elif G.edges[i,j]['bottom'] == 1:
                model.addConstr(gamma[i,j] == 0);
                model.addConstr(gamma[j,t] == 0);
        
    elif restriction == 2:
        for i,j in G.edges:
            if i != s:
                model.addConstr(gamma[i,j] == 0);
    
    
    model.addConstr(gp.quicksum(G.edges[i,j]['cost']*gamma[i,j] for i,j in G.edges) <= budget);
    
    model.addConstr(alpha[t] - alpha[s] >= 1);
    
    for i,j in G.edges:
        model.addConstr(alpha[i] - alpha[j] + theta[i,j] >= 0);
        model.addConstr(beta[i,j] + gamma[i,j] - theta[i,j] >= 0);
        model.addConstr(x[i,j] - G.edges[i,j]['capacity']*(1-gamma[i,j]) <= 0);
            
    model.addConstrs(gp.quicksum(x[v,u] for u in G.successors(v)) -
                     gp.quicksum(x[u,v] for u in G.predecessors(v)) == 0 for v in G.nodes);
       
    
    # ~~~~~~~~~~~~~ Constraints on Number of Victims ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    z = model.addVars(G.edges, vtype=GRB.BINARY, name="z"); 
       
    x1 = model.addVars(G.edges, vtype=GRB.CONTINUOUS, name="x", lb = 0, ub = GRB.INFINITY);
    theta1 = model.addVars(G.edges, vtype=GRB.CONTINUOUS, name="theta", lb = 0, ub = 1); 
    alpha1 = model.addVars(G.nodes, vtype=GRB.CONTINUOUS, name="alpha", lb = 0, ub = 1); 
    beta1 = model.addVars(G.edges, vtype=GRB.CONTINUOUS, name="beta", lb = 0, ub = GRB.INFINITY); 
    
    model.addConstr(x1[t,s] >= gp.quicksum(G.edges[i,j]['capacity']*beta1[i,j] for i,j in G.edges));
    
    model.addConstr(alpha1[t] - alpha1[s] >= 1);
    
    for i,j in G.edges:
        model.addConstr(alpha1[i] - alpha1[j] + theta1[i,j] >= 0);
        model.addConstr(beta1[i,j] + gamma[i,j] - theta1[i,j] >= 0);
        model.addConstr(x1[i,j] - G.edges[i,j]['capacity']*(1-gamma[i,j]) <= 0);
            
    model.addConstrs(gp.quicksum(x1[v,u] for u in G.successors(v)) -
                     gp.quicksum(x1[u,v] for u in G.predecessors(v)) == 0 for v in G.nodes);
    
    for i,j in G.edges:
        if G.edges[i,j]['special'] == 1:
            model.addConstr(x1[i,j] - G.edges[i,j]['capacity']*z[i,j] <= 0);
        else:
            model.addConstr(z[i,j] == 0);
    
    model.addConstr(gp.quicksum(z[i,j] for i,j in G.edges) <= vic);
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    model.addConstr(objVar == x[t,s]);
    
    model.setObjective(gp.quicksum(G.edges[i,j]['capacity']*beta[i,j] for i,j in G.edges), GRB.MINIMIZE);
    

    for i,j in G.edges:
        gamma[i,j].start = Gamma2[i,j];
    
    model.update();
    
    T_Limit = 7200;
    
    model.setParam("IntegralityFocus",1);
    model.setParam('TimeLimit', T_Limit); 
    model.setParam("LogToConsole", 0)
    model.setParam("OutputFlag", 0);
    model.update();
    model.optimize();
    
    Gamma = {};
    
    for i,j in G.edges:
        Gamma[i,j] = gamma[i,j].x;
        
    
    return model.objVal, Gamma;