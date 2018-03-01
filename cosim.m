% FUNCTION FOR CALCULATING COSINE VECTOR SIMILARITY

function sim = cosim(A,B)

sim = dot(A,B)/(norm(A,2)*norm(B,2));