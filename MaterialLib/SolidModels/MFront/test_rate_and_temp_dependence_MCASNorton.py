import mtest

import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['lines.linewidth']= 2.0
plt.rcParams['lines.color']= 'black'
plt.rcParams['legend.frameon']=False
plt.rcParams['font.family'] = 'serif'
plt.rcParams['legend.fontsize']=8
plt.rcParams['font.size'] = 12
# For ipython notebook display set default values.
plt.rcParams['lines.markersize'] = 10
plt.rcParams['grid.linewidth'] = 1
# General settings used by display and print contexts.
#plt.rcParams['axes.axisbelow'] = True
grid_line_color = '0.5'
plt.rcParams['grid.color'] = grid_line_color
plt.rcParams['grid.linestyle'] = '-'


fig, ax = plt.subplots(nrows=2)


for tload in [2,5]:
    for temp in [293,393]:
        sol_t = [0.]
        sol_sxx = [0.]
        sol_syy = [0.]
        sol_p = [0.]
        sol_c = [0.]
        # for theta in [-1.3010636242139548]:
        em = 5.e-2
        steps = [200,100]
        theta = 0
        tmax = 10
        c = np.cos(theta)
        s = np.sin(theta)
        m = mtest.MTest()
        mtest.setVerboseMode(mtest.VerboseLevel.VERBOSE_QUIET)
        m.setMaximumNumberOfSubSteps(10)
        m.setBehaviour('generic', 'src/libBehaviour.so', 'MohrCoulombAbboSloanNorton')
        m.setExternalStateVariable("Temperature", temp)
        m.setImposedStrain('EXX', {0: 0, tload: em * c, tmax: em * c})
        m.setImposedStrain('EYY', {0: 0, tload: em * s, tmax: em * s})
        m.setNonLinearConstraint('SXX+SYY+SZZ', 'Stress')
        m.setNonLinearConstraint('SXY', 'Stress')
        m.setNonLinearConstraint('SXZ', 'Stress')
        m.setNonLinearConstraint('SYZ', 'Stress')
        m.setMaterialProperty('YoungModulus', 150.e3)
        m.setMaterialProperty('PoissonRatio', 0.3)
        m.setMaterialProperty('Cohesion', 3.e1)
        m.setMaterialProperty('FrictionAngle', 30.)
        m.setMaterialProperty('DilatancyAngle', 10.)
        m.setMaterialProperty('TransitionAngle', 29.)
        m.setMaterialProperty('TensionCutOffParameter', 1.e1)
        m.setMaterialProperty('CreepRateFactor', 0.18)
        m.setMaterialProperty('CreepRateExponent', 5.0)
        m.setMaterialProperty('ActivationEnergy', 54.0)
        m.setParameter('UniversalGasConstant', 8.314e-3)
        m.setParameter('ReferenceStress', 1.)
        s = mtest.MTestCurrentState()
        wk = mtest.MTestWorkSpace()
        m.completeInitialisation()
        m.initializeCurrentState(s)
        m.initializeWorkSpace(wk)
        ltime = np.linspace(0,tload,steps[0])
        ltime = np.append(ltime,np.linspace(tload+tload/steps[0],tmax,steps[1]))
        for i in range(len(ltime)-1):
            m.execute(s, wk, ltime[i], ltime[i + 1])
            p = s.getInternalStateVariableValue('EquivalentPlasticStrain')
            c = s.getInternalStateVariableValue('EquivalentCreepStrain')
            sol_t.append(ltime[i+1])
            sol_sxx.append(s.s1[0])
            sol_syy.append(s.s1[1])
            sol_p.append(p)
            sol_c.append(c)
        print("plotting...")
        ax[0].plot(sol_t,sol_sxx,label='$\\sigma_{xx}$, $T = %i$, $t_\\mathrm{load} = %i$' %(temp, tload))
        #ax[0].plot(sol_t,sol_syy,label='$\\sigma_{yy}$, $T = %i$, $t_\\mathrm{load} = %i$' %(temp, tload),ls='--')
        ax[1].plot(sol_t,sol_p,label='$\\epsilon_{pl}$, $T = %i$, $t_\\mathrm{load} = %i$' %(temp, tload))
        ax[1].plot(sol_t,sol_c,label='$\\epsilon_{cr}$, $T = %i$, $t_\\mathrm{load} = %i$' %(temp, tload),ls='--')

print("finishing plot...")
ax[0].set_xlabel('$t$ / d')
ax[0].set_ylabel('$\\sigma$ / MPa')
ax[1].set_xlabel('$t$ / d')
ax[1].set_ylabel('$\\epsilon$')
for i in range(2):
    ax[i].spines['right'].set_visible(False)
    ax[i].spines['top'].set_visible(False)
    ax[i].legend()
    ax[i].set_xlim(0,tmax+5)

fig.tight_layout()
fig.savefig('rate_and_temperature_dependence.pdf')
