def solution(c_0, c_L, D_A, D_B, L_A, L_B, k_plus, k_minus):

    R_A = L_A / D_A
    R_B = L_B / D_B

    phi = (k_plus * c_0 - k_minus * c_L) / (1 + k_plus * R_A + k_minus * R_B)

    c_a_int = c_0 - phi * R_A
    c_b_int = c_L + phi * R_B
    return c_a_int, c_b_int, phi
