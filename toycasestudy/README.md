
# Toy case study instance

The instance of our program represents an interaction graph (IG) and a list of discrete observation. This can be found in  in **instancetoynetwork.lp**.

The predicats used are:

<ul>
<li> <i>sign</i> of arity one which represents the sign that can be taken by a node: {1,-1,0}; 1 is for "+" (overexpressed nodes) -1 for "-" (underexpresssed nodes) and 0 for "0" (unchanged)</li>
<li> <i>edgeS</i> of arity one which represents the sign that can be taken by an edge: {1, -1}; 1 is for "+" (activation) and -1 is for "-" (inhibition)</li>
<li> <i>name</i> of arity one ranging from 1 to k represents the number of artifical influences added in case of inconsistencies, k can be parametrize by the user </li>
<li><i>weight</i> of arity one represents the weight that can be taken by a node: [0..100]</li>
<li><i>observedE</i> of arity three represents the edges of the IG as a triplet composed of: parent node, child node, edge sign.</li>
<li><i>observedV</i> of arity three represents the observed nodes in our IG as a triplet composed of: node, sign, weight. All weights are set to 100 but this could be parametrized. </li>

</ul>
