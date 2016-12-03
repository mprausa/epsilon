(* vim: set syntax=mma expandtab shiftwidth=4 tabstop=4: *)

<<EpsilonTools`

Print["checking transformation.."];

m = Get["matrix.m"];
expected = EpsilonRead["out/epsilon",CheckFuchsian->True,CheckEpsilon->True];

t = Get["out/transformation.m"];
t = t/.EpsilonSymRules[t];

tinv = Inverse[t];

Print[If[And@@(PossibleZeroQ/@Flatten[tinv.m.t - tinv.D[t,x] - expected]),"OKAY :)","FAILED :("]];

