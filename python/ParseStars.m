 (* $Id: ParseStars.m,v 1.6 2014-07-29 20:08:25 oshaughn Exp $
*   Compare to AlbertoCatalogManager
 *)

<<MyTools`TypedList`
<< AppliedStatistics`Nonparametric`PointsToDensity`KernelEstimators`
<< AppliedStatistics`Nonparametric`PointsToDensity`WeightedKernelEstimators`
$dirStarData = $dirROSProjs<>"/PSnucleosynthesis/Branch-GalaxyAssemblyAndZDistribution/Communications/2011-08-03-Jillian-StarData/";

$dirStarDataMore = $dirROSProjs<>"/PSnucleosynthesis/Branch-GalaxyAssemblyAndZDistribution/bhmetals/";

mapfunc = ParallelMap;

<<MyUnits`UnitManager`;

timeScale = 1.223 *10^18 Convert[Second/(Giga Year), 1]


jillianCatalogLoad["Small"] := Module[{tmax,dat},
 SetDirectory[$dirStarData];
 dat = {#[[1]], #[[2]], Random[] }&/@ Drop[Import["starmetals.dat", "Table"], 1];
 tmax = Max[dat[[All,-1]]];
 StarGasCatalog @@ {"Type"-> "Jillian", "ExtractionRedshift"->4, "ExtractionTime"->tmax,
    "Data"->dat
 }
]

jillianCatalogLoad["Large"] := Module[{tmax,dat},
 SetDirectory[$dirStarData];
 dat = {#[[1]], #[[2]], #[[3]] timeScale}&/@ Drop[Import["starmetals.20110819.dat", "Table"], 1];
 tmax = Max[dat[[All,3]]];
 StarGasCatalog @@ {"Type"-> "Jillian", "ExtractionRedshift"->4, "ExtractionTime"->tmax,
    "Data"->dat }
];

jillianCatalogLoad["Large","Movie", z_]/;MemberQ[{2, 2.4, 2.6, 3, 3.4, 4,5,6,8},z] := Module[{tmax,dat},
 SetDirectory[$dirStarDataMore];
 dat = {#[[1]], #[[2]], #[[3]] timeScale}&/@ Drop[Import["h258.z" <>ToString[z]<>".starmetal.dat.gz", "Table"], 1];
 tmax = Max[dat[[All,3]]];
 StarGasCatalog @@ {"Type"-> "Jillian", "ExtractionRedshift"->z, "ExtractionTime"->tmax,
    "Data"->dat
 }
]

jillianCatalogLoad["Large","H2CoolingOn"]:= Module[{tmax,dat},
 SetDirectory[$dirStarDataMore];
 dat = {#[[1]], #[[2]], #[[3]] timeScale}&/@ Drop[Import["h285.starmetal.H2.dat.gz", "Table"], 1];
 tmax = Max[dat[[All,3]]];
 StarGasCatalog @@ {"Type"-> "Jillian", "ExtractionRedshift"->z, "ExtractionTime"->tmax,
    "Data"->dat
 }
]

jillianCatalogLoad["Large","H2CoolingOn", "Movie",z_] /;MemberQ[{2,3, 4, 4.8, 5, 6, 8},z] := Module[{tmax,dat},
 SetDirectory[$dirStarDataMore];
  If[z>2,
  (* Ages are in Gyr already *)
    dat =  Drop[Import["h258.H2.z" <>ToString[z]<>".starmetal.dat.gz", "Table"], 1];,
    dat = {#[[1]], #[[2]], #[[3]] timeScale}&/@ Drop[Import["h285.starmetal.H2.dat.gz", "Table"], 1];
 ];
 tmax = Max[dat[[All,3]]];
 StarGasCatalog @@ {"Type"-> "Jillian", "ExtractionRedshift"->z, "ExtractionTime"->tmax,
    "Data"->dat
 }
]




jcReduce[jc_StarGasCatalog, cond_] := Module[{dat},
   dat = Select[jc[[-1,2]], cond];
 StarGasCatalog @@ {"Type"-> "Jillian", "ExtractionRedshift"->4, "ExtractionTime"->"XX",
    "Data"->dat
 }
]

jcEvaluateAndWeight[jc_StarGasCatalog, expr_, exprWeight_, Zg_, Zs_, t_] := Module[{dat},
  dat = jc[[-1,2]];  (* assume type works *)
 
    mapfunc[({expr,  exprWeight} /. { Zs->#[[1]], Zg->#[[2]], t->#[[3]]}
    )& , dat]

]

jcIntegrateExpression[jc_StarGasCatalog, expr_, exprWeight_, Zg_, Zs_, t_] := Module[{},
  Plus @@ (Times[Sequence @@ #] & /@ jcEvaluateAndWeight[jc, expr, exprWeight,Zg, Zs, t])
]



jcFastSFR[cat_StarGasCatalog, opts___] := Module[{mstarParticles},
   mstarParticles = 26676*0.3.;  (* 0.3 of particle mass ... a little drift, but usually ok. *)
  tAgeLimit = "AgeLimitGyr" /. {opts} /. "AgeLimitGyr" -> 0.02;
     datReduced =  Select[cat[[-1, 2]], #[[-1]] < (1.0) (tAgeLimit) &]; 
    Length[datReduced]*mstarParticles/(10^9*tAgeLimit)   
];
jcFastMstar[cat_StarGasCatalog, opts___] := Module[{mstarParticles},
   mstarParticles = 26676*0.3.;  (* 0.3 of particle mass *)
      datReduced =  cat[[-1, 2]];
    Length[datReduced]*mstarParticles
];


(* Weighting by star mass.  These are all the same, so easy *)
jcFastLogZCumulativePlot[cat_StarGasCatalog, opts___]:= Module[
  {tMax, tAgeLimit, datReduced, datRaw, fn, nIndex, tAgeRange},
   tMax = Max[ cat[[-1, 2]][[All, 3]]];
   nIndex = "Index"/.{opts}/."Index"->1;
   tAgeLimit = "AgeLimitGyr" /. {opts} /. "AgeLimitGyr" -> tMax;
   tAgeRange = "AgeRangeGyr"/.{opts}/."AgeRangeGyr"->Null;

    If[ SameQ[tAgeRange, Null] || !MatchQ[tAgeRange, _List],
     datReduced =  Select[cat[[-1, 2]], #[[-1]] < (1.0) (tAgeLimit) &]; ,
     datReduced =  Select[cat[[-1, 2]], #[[-1]] <  (tAgeRange[[2]] ) && #[[-1]] >  (tAgeRange[[1]] ) &]; 
    ];

     datRaw = Estimate1dPcumWithHardStepViaSortRaw[Log[10, datReduced[[All,nIndex]]//Union  ]];
     fn = Interpolation[datRaw, InterpolationOrder->1];
     Plot[Log[10, fn[x]], {x,-7, -1}, Evaluate[Sequence@@FilterRules[{opts}, Options[Plot] ]] ]
]
jcGaussianModelLogZ[cat_StarGasCatalog, lZcut_, opts___] := 
 Module[{tMax, tAgeLimit, datReduced, fn, lPCut}, 
  tMax = Max[cat[[-1, 2]][[All, 3]]];
  nIndex = "Index" /. {opts} /. "Index" -> 1;
  tAgeLimit = "AgeLimitGyr" /. {opts} /. "AgeLimitGyr" -> tMax;
  lPCut = "LogProbabilityCut"/.{opts}/."LogProbabilityCut"-> -3;
  datReduced = 
   Select[cat[[-1, 
      2]], #[[-1]] < (1.01) (tAgeLimit) && #[[nIndex]] > 10^lZcut &];
  datRaw = 
   Drop[Estimate1dPcumWithHardStepViaSortRaw[
     Log[10, datReduced[[All, nIndex]] // Union]], 1];
  datRaw = Select[datRaw, #[[2]] > 10^lPCut &];  (* throw out insignificantly small probability events that aren't gaussian *)
  myfit = 
   NonlinearModelFit[{#[[1]], Log[10, #[[2]]] // N} & /@ datRaw, 
    Log[10, CDF[NormalDistribution[lZ, sZ]][x]], {lZ, sZ}, x];
  {myfit["ParameterTableEntries"][[All, 1]], myfit}
  ]


(* DEPRECATED *)
showHistIndx[cat_, indx_, opts___] := Module[{tAgeLimit, tMax, fMax},
   tMax = Max[ cat[[-1, 2]][[All, 3]]];
   tAgeLimit = "AgeLimitGyr" /. {opts} /. "AgeLimitGyr" -> tMax;
   fMax = 
    "AgeLimitFractional" /. {opts} /. "AgeLimitFractional" -> Null;
   If[! SameQ[fMax , Null] && fMax < 1 && fMax > 0, 
    tAgeLimit = fMax tMax;
    ];
   
   datReduced = 
    Select[cat[[-1, 2]], #[[-1]] < (1.01) (tAgeLimit) &]; 
   Print[Length[datReduced]];
   
   Histogram[If[! SameQ[indx, 3], 
     Log[10,
      datReduced[[All, indx]] // Sort
      ],
     datReduced[[All, indx]]
     ], If[! SameQ[indx, 3], {0.1}, {0.001}], {"Log", "CDF"}, 
    Axes -> True, PlotRange -> All, 
    Sequence @@ FilterRules[{opts}, Options[Histogram]], 
    Frame -> True, FrameLabel -> If[indx < 3,
      {"log Z", "log P(<Z)"},
      {"log t", "dP/dlog t"}
      ]]];
