filepath = NotebookDirectory[]                                      (* files weight spectrum.nb, data.txt и result.txt are in the same folder *)
basis = ToExpression[(Characters@DeleteDuplicates[
        Import[FileNameJoin[{filepath, "data.txt"}], "Lines"]])];
ws = WeightSpectrumGrayCode[basis]

(* compilation of the weight spectrum calculation *)
(* WeightSpectrum=
   Compile[{{basevectors,_Integer,2}},
           Module[{
                          n=Dimensions[basevectors][[-1]],          (* dimension of vectors *)
                          k=Length[basevectors],                    (* number of vectors *)
                          s=Table[0,{i,0,Dimensions[basevectors][[-1]]}]      
                  },
                  Do[
                          With[{l=Total[Mod[Total[IntegerDigits[b,2,k]*basevectors],2]]+1},
                          s[[l]]=s[[l]]+1],
                          {b,0,2^k-1}];
                 Return[Transpose[{Range[0,Last[Dimensions[basis]]],s}]]
           ]
   ];   *)

(* compiling the weight spectrum calculation using Gray code *)
WeightSpectrumGrayCode = Compile[{{basevectors, _Integer, 2}},
   Module[{
       n = Dimensions[basevectors][[-1]],                           (* dimension of vectors *)
       k = Length[basevectors],                                     (* number of vectors *)
       s = Table[0, {n + 1}],   
       currentVector = Table[0, {n}],
       m = 0, l = 0                                                 (* the first combination is always {0,0,..} *)
   },
   
       s[[1]] = 1;
       Do[
           m = Log2[BitAnd[-1 - b, b + 1]] + 1;                     (* changed bit position *) 
           currentVector = BitXor[currentVector, basevectors[[m]]]; (* there is no complete addition of vectors, but only the addition of the difference *) 
           l = Total[currentVector] + 1;                            (* weight count *)
           s[[l]] = s[[l]] + 1,
           {b, 0, 2^k - 2}
       ];
           Return[Transpose[{Range[0, Last[Dimensions[basis]]], s}]]
   ]
 ]; 
   
(* we divide the range of all combinations into equal parts depending on the number of available cores,
that is, as a result, each kernel returns a partial spectrum *)
partition[{i1_Integer, iN_Integer}, n_Integer] :=
    With[{dn = Round[(iN - i1 + 1)/n]},
        Join[
            {{i1, i1 + dn - 1}},
            Table[{dn*(i - 1), i*dn - 1}, {i, 2, n - 1}],
            {{(n - 1)*dn, iN}}
        ]
    ]
  
(* send these calculations to different cores *)
WeightSpectrumCodeParallel[basis : {{__Integer} ..}] :=
    With[{
        kernels = (If[Kernels[] === {}, LaunchKernels[]]; Kernels[]),
        parts = partition[{0, 2^Length[basis] - 1}, Length[Kernels[]]]
    },
        Total[WaitAll[Table[
             ParallelEvaluate[ParallelSubmit[WeightSpectrumGrayCode[basis, parts[[$KernelID]]]], kernel],
             {kernel, kernels}
        ]]]
    ]
  
(* render computation time *)
AbsoluteTiming[ParallelEvaluate[WeightSpectrumGrayCode[basis]];][[1]]

(* show all running kernels *)
Kernels[]

(* save result to file result.txt *)
Export[
       FileNameJoin[{filepath, "result.txt"}], ws, "Table", "FieldSeparators" -> "\t"]

(* look at the result in the file *)
SystemOpen[%]

(* visualize the weight spectrum *)
ListPlot[ws, PlotTheme -> "DarkColor", AxesLabel -> {"Вес", "Количество"}, AxesOrigin -> {0, 0}, PlotRange -> All, Filling -> Bottom]
