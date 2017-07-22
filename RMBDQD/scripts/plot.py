#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys

nfile = sys.argv[1]

E, w, dp = np.genfromtxt(nfile, unpack=True)

plt.figure()

plt.plot(E,[0 for i in E])
plt.plot(E,dp)
# plt.title(u'$100$ pontos')
if max(E)<0.5:
    plt.xlabel(u'Coeficiente de Poisson')
else:
    plt.xlabel(u'Módulo de Elasticidade (Pa)')
plt.ylabel(u'Desvio Padrão $(\omega_{ exp } - \omega_{fem})$')

plt.xlim(min(E), max(E))
# plt.plot(E,w)

plt.grid()
plt.show()
