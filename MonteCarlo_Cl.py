import numpy as np
import powerspec_HOD as hod


N = 1000
#limit

M_halo = 10.**(10 + (8.*np.random.rand(N)))
z4inte = 0.2 + (0.4*np.random.rand(N))

N_cent, N_sat= hod.N_hod_3param(M_halo, 1e11, 1e13, 0.7)
print(1/N *sum(hod.__dndm_func(M_halo, 0.2)*(N_cent+N_sat)))

