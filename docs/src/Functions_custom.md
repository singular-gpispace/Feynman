```@meta
CurrentModule = Feynman
```


# FUNCTIONS TO WORK WITH USER DEFINED PROPAGATORS

```@raw html
<details>
<summary>makeBaikovMatrix(def G,list internalmomenta,list externalmomenta,list mandelsonvars,list propagat,list replacementRules)</summary>
```
**USAGE**   :  makeBaikovMatrix(G,internalmomenta,externalmomenta, mandelsonvars, propagat, replacementRules); G labeledgraph, or G graph

**ASSUME**  :   G is a Graph, or@*
                G is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.

**RETURN**  :  a labeled graph G1, computes the Baikov matrix of G defined in G1.baikovover and stores it in G1.baikovmatrix

**Example** :

```singular
//Setting the graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  //include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,t1,t2),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(t1,t2);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*t1,p1*p4,1/2*t2,p2*p4,-1/2*(t1+t2));
  
  //compute Baikov matrix
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: makeM1(def G0)</summary>
```
**USAGE**   :  makeM1(G0); G labeledgraph, or G graph

**ASSUME**  :   0 is a Graph, or@*
                G0 is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.

**RETURN**  :  The module M1 over G1.baikovover that requires to compute IBP identities 

**Example** :

```singular
    graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
    labeledgraph G1=computeBaikovMatrix(G);
    ring RB=G1.baikovover;
    RB;
    module ML=makeM1(G1);
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: makeM2(def G0,list Nu)</summary>
```
**USAGE**   :  pcomputeM2(G,Nu); G labeledgraph, or G graph

**ASSUME**  :   G is a Graph, or@*
                G is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.
                Nu is the seed.

**RETURN**  :  The module M2 over G1.baikovover that requires to compute IBP identities

**Example** :

```singular
    graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
    labeledgraph G1=computeBaikovMatrix(G);
    ring RB=G1.baikovover;
    RB;
    module M2=makeM2(G1,list(1,1,1,0,0,1,0,0,0));
    M2;
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: makeFormalIBP(def G0,list sector)</summary>
```
**USAGE**   :  makeFormalIBP(G0,sector); G0 graph

**ASSUME**  :   G0 is the labelled graph and is the output of makeBaikovMatrix. 
                sector is the list of integers represent a sector of G.

**RETURN**  :  generators of the standard basis of  the module M1 intersect M2.

**Example** :

```singular
//graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  //include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,t1,t2),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(t1,t2);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*t1,p1*p4,1/2*t2,p2*p4,-1/2*(t1+t2));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  ring RZ=G1.baikovover;
  list sector=list(1,2,3); //sector that we are interested
  module M=makeFormalIBP(G1,sector);
  size(M);
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: makeIBPVec(def G0,def M12, list setNu)</summary>
```
**USAGE**   :   makeIBPVec(G0,M12,setNu); G0 graph, M12 module, setNu list,

**ASSUME**  :   setNu is a list of seed correspond to the graph G0 which are belong to the same sector 
                M12 is the formal IBP of the corresponding sector (output of makeFormalIBP).

**RETURN**  :  setIBP S, where it contains all the IBP relations obtained by module intersection and seeding

**Example** :

```singular
 //include graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  //include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,t1,t2),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(t1,t2);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*t1,p1*p4,1/2*t2,p2*p4,-1/2*(t1+t2));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  ring RZ=G1.baikovover;

  //Assume the seed belong to the same sector
  list sector=list(1,2,3,4,5,6,7);
  list seeds=list(list(1,1,1,1,1,1,1,0,0));

  module M=makeFormalIBP(G1,sector);
  setIBP S=makeIBPVec(G1,M,seeds);

    oneIBP I=S.IBP[1];
    I;
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: getSortMeasuresVec(vector l,int x)</summary>
```
**USAGE**   :  getSortMeasures(l,x); l list, x int; 

**ASSUME**  :  l is a list of integers (i.e a seed) and x is number of Baikov variables.

**RETURN**  :  a vector of sort measures that are used in Laporta Algorithm

```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: getSortedIntegralsVec(setIBP I)</summary>
```
**USAGE**   :  etSortedIntegrals(I); I setIBP

**ASSUME**  :  

**RETURN**  :   list ind where each entry is a pair (indv,sortmeasures),
                indv is the list of indices(seed) appered in the setIBP 
                and sortmeasures is the output of getSortMeasuresVec(indv,x).
                The function getSortedIntegrals extract the seeds appeared in the IBP identities of the setIBP,
                sort them lexicographically based on the values got from getSortMeasuresVec and return the output.

**Example** :

```singular
//include graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  //include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,t1,t2),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(t1,t2);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*t1,p1*p4,1/2*t2,p2*p4,-1/2*(t1+t2));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  
  ring RZ=G1.baikovover;
  list sector=list(1,2,3,4,5,6,7);
  list seeds=list(list(1,1,1,1,1,1,1,-4,0),list(1,1,1,1,1,1,1,-1,-3),list(1,1,1,1,1,1,1,-2,-2),list(1,1,1,1,1,1,1,-3,-1),list(1,1,1,1,1,1,1,0,-4));
 
  module M=makeFormalIBP(G1,sector);
  setIBP S=makeIBPVec(G1,M,seeds);
  list L =getSortedIntegralsVec(S);
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: extractCoefVec(oneIBP I,list ind,list l,list sector)(matrix M)</summary>
```
**USAGE**   :  extractCoefVec(I,ind,l,sector); I oneIBP,ind list,l list

**ASSUME**  :   ind is the output of getSortedIntegralsVec, and l is the list of values over the base field I.baikovover. 
                size(l)=npars(I.baikovover)

**RETURN**  :   list of values where, the i-th element is the evaluation of coefficient function  at values in the list l 
                of the IBP relation oneIBP, whose index is i=ind[i][1].

```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: makeMatVec(setIBP S,list val,list ind,list sector)</summary>
```
**USAGE**   :  makeMatVec(S,val,ind,sector); S setIBP,ind list,l list,sector list;

**ASSUME**  :  size(val)=npars(S.over), ind is the output of getSortedIntegralsVec

**RETURN**  :   matrix,where i-th row correspond to the evaluation of coefficient functions of i-th IBP in setIBP. 
                Columns of the matrix correspond to the all used indices in the setIBP which are ordered with 
                respect to the output ofgetSortMeasures.

**Example** :

```singular
// include graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  //include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,t1,t2),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(t1,t2);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*t1,p1*p4,1/2*t2,p2*p4,-1/2*(t1+t2));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  
  ring RZ=G1.baikovover;

  list sector=list(1,2,3,4,5,6,7);
  list seeds=list(list(1,1,1,1,1,1,1,-4,0),list(1,1,1,1,1,1,1,-1,-3),list(1,1,1,1,1,1,1,-2,-2),list(1,1,1,1,1,1,1,-3,-1),list(1,1,1,1,1,1,1,0,-4));
  
  list val=getRandom(93187,npars(RZ));
  
  module M=makeFormalIBP(G1,sector);
  setIBP S=makeIBPVec(G1,M,seeds);
  
  list L =getSortedIntegralsVec(S);
  matrix N=makeMatVec(S,val,L,sector);
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: getReducedIBPVec(setIBP S,int p,list sector)</summary>
```
**USAGE**   :  getRedIBPVec(S,p,sector); 

**ASSUME**  :  S is setIBP, and p is a prime number. 

**RETURN**  :  list L, L[1]=indIBP, L[2]=seed where,
                indIBP contain the linearly independent IBP relations of setIBP which are obtained by finite 
                field row reduction over the field Fp. 
                seed contain the indeces correspond to the non-free columns in rref.

**Example** :

```singular
graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  ring R=0,(q1,q2,p1,p2,p4,s,t),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(s,t);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*s,p1*p4,1/2*t,p2*p4,-1/2*(s+t));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  ring RZ=G1.baikovover;
  list sector=list(1,2,3,4,5,6,7);
  list seeds=list(list(1,1,1,1,1,1,1,-4,0),list(1,1,1,1,1,1,1,-1,-3),list(1,1,1,1,1,1,1,-2,-2),list(1,1,1,1,1,1,1,-3,-1),list(1,1,1,1,1,1,1,0,-4));
  module M=makeFormalIBP(G1,sector);
  setIBP S=makeIBPVec(G1,M,seeds);

  list L=getReducedIBPVec(S,93187,sector);
  size(L[1])<size(S.IBP);
  ring RS=S.over;
  MI mi=L[2];
  size(mi.masterIntegrals);
  mi;
  //print of reduced IBPs
  for(int i=1;i<=size(L[1]);i++){
    oneIBP I=L[1][i];
    I;
  }
  //print of all IBPs
  for(int i=1;i<=1;i++){
    oneIBP I=S.IBP[i];
    I;
  }
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: getReleventIBPs(setIBP S,def sector)</summary>
```
**USAGE**   :  getReleventIBPs(S,sector); 

**ASSUME**  :  S is setIBP. 

**RETURN**  :  setIBP S, where for each IBP, the terms consist of the integrals that are not 
                belong to the given sector are removed (i.e., masking process is imposed).

**Example** :

```singular
//include graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

//include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,s,t),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(s,t);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*s,p1*p4,1/2*t,p2*p4,-1/2*(s+t));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  ring RZ=G1.baikovover;
  list sector=list(1,2,3,4,5,6,7);
  list seeds=list(list(1,1,1,1,1,1,1,-4,0),list(1,1,1,1,1,1,1,-1,-3),list(1,1,1,1,1,1,1,-2,-2),list(1,1,1,1,1,1,1,-3,-1),list(1,1,1,1,1,1,1,0,-4));
  
  
  module M=makeFormalIBP(G1,sector);
  
  setIBP S=makeIBPVec(G1,M,seeds);
  setIBP S1=getReleventIBPs(S,sector);

  oneIBP I = S1.IBP[1];
  I;
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: getProperIBPs(int nv, list L)</summary>
```
**USAGE**   :  getProperIBPs(nv,L); 

**ASSUME**  :  L is a list of IBPs and nv is number of Baikov variables. 

**RETURN**  :  a list of IBPs where the indeces of integrals (in vector format) in IBPs are converted to lists.

**Example** :

```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: getReducedIBPwithMask(def G1,def M,int p,list sector)</summary>
```
**USAGE**   :  getReducedIBPwithMask(G1,M,p,sector); 

**ASSUME**  :  G is a graph, M is the output of the function makeFormalIBP, p is a prime number and sector 
                is a list. 

**RETURN**  :   a list where the entry i contain the list of independent IBPs correspond to step i-1. Here we 
                consider steps upto 4.

**Example** :

```singular
    graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
    labeledgraph lG = labelGraph(G,0);
    labeledgraph G1 = eliminateVariables(lG);
    labeledgraph G2 = removeElimVars(G1);

    ring R=0,(q1,q2,p1,p2,p4,s,t),dp;
    list internalmomenta=list(q1,q2);
    list externalmomenta=list(p1,p2,p4);
    list mandelsonvars=list(s,t);
    list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
    list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*s,p1*p4,1/2*t,p2*p4,-1/2*(s+t));
  
    labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
    ring RZ=G1.baikovover;
    list sector=list(1,2,3,4,5,6,7);
    module M=makeFormalIBP(G1,sector);
    list  Lprop=getReducedIBPwithMask(G1,M,93187,sector);  
    //One can print the independent IBPs with mask correspond to step 0.
    for(int j=1;j<=size(Lprop[1]);j++){
     oneIBP I=Lprop[1][j];
        if(size(I.i)<>0){
        print(I);
        }
    } 
```
```@raw html
</details>
```
```@raw html
<details>
<summary>procedure: getIBPwithMask(def G1,def M,int p,list sector)</summary>
```
**USAGE**   :  getIBPwithMask(G1,M,p,sector);

**ASSUME**  :  G is a graph, M is the output of the function makeFormalIBP, p is a prime number and sector
                is a list.

**RETURN**  :   a list where the entry i contain the list of IBPs correspond to step i-1. Here we consider 
                steps upto 4.

**Example** :

```singular
//include graph information
  graph  G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
  labeledgraph lG = labelGraph(G,0);
  labeledgraph G1 = eliminateVariables(lG);
  labeledgraph G2 = removeElimVars(G1);

  ////include user specified propagators, replacement rules etc.
  ring R=0,(q1,q2,p1,p2,p4,s,t),dp;
  list internalmomenta=list(q1,q2);
  list externalmomenta=list(p1,p2,p4);
  list mandelsonvars=list(s,t);
  list propagat=list(q1^2,(q1-p1)^2,(q1-p1-p2)^2,(q2+p1+p2)^2,(q2-p4)^2,q2^2,(q1+q2)^2,(q1+p4)^2,(q2+p1)^2);
  list replacementRules=list(p1^2,0,p2^2,0,p4^2,0,p1*p2,1/2*s,p1*p4,1/2*t,p2*p4,-1/2*(s+t));
  
  labeledgraph G1=makeBaikovMatrix(G2,internalmomenta,externalmomenta,mandelsonvars,propagat,replacementRules);
  ring RZ=G1.baikovover;
  list sector=list(1,2,3,4,5,6,7);
  module M=makeFormalIBP(G1,sector);
  list Lprop=getIBPwithMask(G1,M,93187,sector);  

  //One can use the following to print the IBPs with mask correspond to step 0
   for(int j=1;j<=size(Lprop[1]);j++){
    oneIBP I=Lprop[1][j];
   
      print(I);
    
   } 
```
```@raw html
</details>
```
