function dif = Difference(k, lon_stiff, ben_stiff)
global ForceMat
global BeamMat
BeamMat = BeamForce(k, lon_stiff, ben_stiff);
dif = ForceMat - BeamMat;
end