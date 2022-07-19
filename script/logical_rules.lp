
%-------------Author: Sophie Le Bars-------------------
%-------Allow to predict unobserved nodes in a network and give out weight associated with this prediction
% input: a network(nodes: vertex, edges: observedE) with some observed nodes inside our network (observedV) 
% output: unobserved nodes are predicted in the graph


sign(1;-1;0).%All possible signs for observed or predicted nodes.
edgeS(1;-1).
name(1..k).

weight(0..100). %All possible weights for observed or predicted nodes.





%%%%%%%%%%%%%%%%%%%%Possible complete weighted labellings (Section 3.2 le bars et al. BMC supplements)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


1{labelV(I,S,W) : sign(S) , weight(W)}1 :- vertex(I).%Generate a weighted label (labelV) for each node in the IG.


%%%%%%%%%%%%%%%%%%%%Management of input nodes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


input(P):-vertex(P),not observedE(_,P,_).%Retrieve nodes of the network without parent nodes: input node.

:- input(I),labelV(I,S,W), W!=100.%Fix the weight of the input node to 100.


%%%%%%%%%%%%%%%%%%%%Test rules MajoritySign (Section 3.1.2 le bars et al. BMC supplements)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%Experimental observation signs are kept%%%


observedV2(I,S) :- observedV(I,S,_).

1{labelV(I,S,W) : sign(S) , weight(W)}1 :- observedV2(I,S).%The observed nodes keep their sign so their weighted label too


%%%Signed Majority wins%%%


parent(P,E,100):- observedE(P,E,Y),labelV(P,S,X).%Get for a node its parents' nodes (parent: parent node, children node, weight). The weight is fixed to 100 for the propagation and could be replaced by X, the parent's weight. 

parentsigne(P,E,S*Y):- observedE(P,E,Y),labelV(P,S,X).%Get for a node its parents' nodes received influences (parentsigne: parent node, children node, received influence). The influence is defined as a multiplication of the parent's sign and the edge's sign.  

countsign(I,X,Y):-vertex(I),sign(X),Y=#count{P:parentsigne(P,I,X)},not input(I).%Count for a node the number of parents nodes giving the same influences between three possible received influences (-1, +1, 0). The node will have a number associated (Y) for each possible received influence. If it doesn't receive any influence, this number will be 0. (countsign: node, received influence sign, number of received influences for the sign).

signmaj(I,X) :- vertex(I), countsign(I,X,S1), countsign(I,Y,S2), S1>S2,X!=0,Y!=0.%For a node (I) gives its majoritarian sign between -1 and 1 (0 are managed in Balanced).


%%%Balanced%%%


labelV(I,0,100) :-vertex(I); not signmaj(I,_), countsign(I,1,Z),countsign(I,-1,Z), countsign(I,0,X),not forbidden(I,0), not fixed(I).%For a node, if the number of received "+" influence is equal to the number of received "-" influence (countsign, Z), then it will be signed as "0".

labelV(I,0,100) :- vertex(I); not signmaj(I,_), countsign(I,0,Z),Z=#count{P:parent(P,I,_)} ,not forbidden(I,0), not fixed(I).%For a node, if all the received influences are "0" influences, it will be signed as "0".


%%%Weight assignment%%%

sumweight(I,X,Y):-vertex(I),sign(X),Y=#sum{Z,P:parent(P,I,Z),parentsigne(P,I,X)},not input(I).%Sum for a node the total weight associated with its received influences (-1,+1,0). (sumweight: node, received influence sign, the total weight associated with the received influence sign).


labelV(I,S,W) :- signmaj(I,S); sumweight(I,S,N) ; C=#count{P:parent(P,I,_)};W=N/C; S!=0; not forbidden(I,S);not fixed(I).%If the node is consistent and not repaired (not forbidden, not fixed) and has a majoritarian sign (S:{1,-1}). Its weight is equal to the ratio between the sum of the parent's weights, holding the majoritarian sign, and its number of parents. 
labelV(I,S,0) :- labelV(I,S,_); signmaj(I,S); fixed(I);S!=0; not forbidden(I,S).%If a node with a majoritarian sign is repaired, then its weight is fixed to 0.

:- labelV(I,S,0),not signmaj(I,S),S!=0.%We cannot have a weighted label with a sign not equal to 0 if it does not have a majoritarian sign. 
:- labelV(I,S,W);not signmaj(I,S); fixed(I); weight(W);W!=0.%If a node with a sign equal to  0 is repaired, its weight is fixed to 0.



%%%%%%%%%%%%%%%%%%%%Consistency (Section 3.1.3 le bars et al. BMC supplements)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


forbidden(N,S) :- vertex(N); not input(N); not signmaj(N,S);observedV2(N,S); S!=0.%If a node is observed, its majority sign cannot be different to its observed sign, so we forbid this event to occur with the predicate forbidden. (forbidden: node, sign).

nodeforbidden(I,S) :- labelV(I,S,_); forbidden(I,S).%A node is forbidden if the sign of its weighted label is forbidden.

pred(I,S,W):- labelV(I,S,W); not forbidden(I,S).%A node is consistent if its weighted label is not forbidden.


:- nodeforbidden(X,S),sign(S),vertex(X).%This line prevents node forbidden from occurring in the proposed answer sets.


%%%%%%%%%%%%%%%%%%%%Repairs (Section 3.1.3 le bars et al. BMC supplements)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


new_influence(X,1) :- not not new_influence(X,1), not new_influence(X,-1), vertex(X).%Generate new negative influence in the graph to reestablish consistency. Same line as Iggy.
new_influence(X,-1) :- not not new_influence(X,-1), not new_influence(X,1), vertex(X).%Generate new positive influence in the graph to reestablish consistency. Same line as Iggy.



:- new_influence(X,_),input(X).%We don't want the input node to have new influence over it.

:- new_influence(X,S), forbidden(X,T), S!=T,sign(T).%Only generate new influence, which explains the sign of a forbidden node.

parentart(P,E):-parentsigne(P,E,_),name(P).%Predicate used to minimise new influence added. (parentart: artificial parent, children node).


fixed(I) :- vertex(I), parentsigne(P,I,_),name(P).%Predicate to know if a node is fixed. (fixed: node).




0{parentsigne(P,E,S):name(P)}k:- new_influence(E,S), vertex(E).%Add artificial nodes to the network from 0 to k per node present in the network.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimisation (Section 3.1.4 le bars et al. BMC supplements)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#minimize{
           1,(X,E)  : parentart(X,E)%Minimize added artificial influences per node. 

         }.


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Print results (optimal answer sets Section 3.1.4 le bars et al. BMC supplements)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#show pred/3. %Weighted label of nodes in the graph

%#show fixed/1. %Repaired nodes (uncomment this line to get them)


