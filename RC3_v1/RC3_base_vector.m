function x = RC3_base_vector(P)
%RC3_BASE_VECTOR  Decision vector of the base case (built by RC3_setup):
%tie switches open, PCC closed (unless forced open by the island mode), no
%DG / SC / transformer / D-STATCOM / UPFC / EVCS / BSS installed, OLTC at
%tap 0, no demand response.
x = P.xbase;
end
