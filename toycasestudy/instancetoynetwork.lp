%%%%%%%%%% Instance of toy network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sign(1;-1;0).%All possible signs for observed or predicted nodes
edgeS(1;-1). %Sign for edges 1 = activation; -1 = inhibition

name(1..k). %The number of artificial influences added k is parametrized in our toy example k=2

weight(0..100).%All possible weights for observed or predicted nodes

%%%%%%%%% Interaction Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex(i;b;e;c;g;a;h;d;f;j).%The 10 nodes of the IG toy case study


%The observed edges of the IG toy case study: first the parent node, second the child and third the sign of the edges.
observedE(e,d,-1).
observedE(f,g,1).
observedE(d,f,1).
observedE(c,d,1).
observedE(b,d,1).
observedE(h,g,1).
observedE(a,d,1).
observedE(i,g,1).
observedE(j,i,1).


%%%%%%%%%% Experimental observation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The observed nodes: first the node name, second its sign and third its associated weight
observedV(f,1,100).
observedV(b,1,100).
observedV(a,-1,100).
observedV(c,1,100).
observedV(h,0,100).
observedV(i,1,100).
observedV(j,-1,100).




