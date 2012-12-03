import numpy as np
import __builtin__
pymin = __builtin__.min
pymax = __builtin__.max

class LineSearch:
    def __init__(self,  xtol=1e-14):

        self.xtol = xtol
        self.task = 'START'
        self.isave = np.zeros((2,), np.intc)
        self.dsave = np.zeros((13,), float)
        self.fc = 0
        self.gc = 0
        self.case = 0
        self.old_stp = 0

    def _line_search(self, func, myfprime, xk, pk, gfk, old_fval, old_old_fval,
                     maxstep=.2, c1=.23, c2=0.46, xtrapl=1.1, xtrapu=4.,
                     stpmax=50., stpmin=1e-8, args=()):
        self.stpmin = stpmin
        self.pk = pk
        p_size = np.sqrt((pk **2).sum())
        self.stpmax = stpmax
        self.xtrapl = xtrapl
        self.xtrapu = xtrapu
        self.maxstep = maxstep
        phi0 = old_fval
        derphi0 = np.dot(gfk,pk)
        self.dim = len(pk)
        self.gms = np.sqrt(self.dim) * maxstep
        #alpha1 = pymin(maxstep,1.01*2*(phi0-old_old_fval)/derphi0)
        alpha1 = 1.
        self.no_update = False

        if isinstance(myfprime,type(())):
            eps = myfprime[1]
            fprime = myfprime[0]
            newargs = (f,eps) + args
            gradient = False
        else:
            fprime = myfprime
            newargs = args
            gradient = True

        fval = old_fval
        gval = gfk
        self.steps=[]

        while 1:
            stp = self.step(alpha1, phi0, derphi0, c1, c2,
                                             self.xtol,
                                             self.isave, self.dsave)

            if self.task[:2] == 'FG':
                alpha1 = stp
                fval = func(xk + stp * pk, *args)
                self.fc += 1
                gval = fprime(xk + stp * pk, *newargs)
                if gradient: self.gc += 1
                else: self.fc += len(xk) + 1
                phi0 = fval
                derphi0 = np.dot(gval,pk)
                self.old_stp = alpha1
                if self.no_update == True:
                    break
            else:
                break

        if self.task[:5] == 'ERROR' or self.task[1:4] == 'WARN':
            stp = None  # failed
        return stp, fval, old_fval, self.no_update

    def step(self, stp, f, g, c1, c2, xtol, isave, dsave):
        if self.task[:5] == 'START':
            # Check the input arguments for errors.
            if stp < self.stpmin:
                self.task = 'ERROR: STP .LT. minstep'
            if stp > self.stpmax:
                self.task = 'ERROR: STP .GT. maxstep'
            if g >= 0:
                self.task = 'ERROR: INITIAL G >= 0'
            if c1 < 0:
                self.task = 'ERROR: c1 .LT. 0'
            if c2 < 0:
                self.task = 'ERROR: c2 .LT. 0'
            if xtol < 0:
                self.task = 'ERROR: XTOL .LT. 0'
            if self.stpmin < 0:
                self.task = 'ERROR: minstep .LT. 0'
            if self.stpmax < self.stpmin:
                self.task = 'ERROR: maxstep .LT. minstep'
            if self.task[:5] == 'ERROR':
                return stp

            # Initialize local variables.
            self.bracket = False
            stage = 1
            finit = f
            ginit = g
            gtest = c1 * ginit
            width = self.stpmax - self.stpmin
            width1 = width / .5
#           The variables stx, fx, gx contain the values of the step,
#           function, and derivative at the best step.
#           The variables sty, fy, gy contain the values of the step,
#           function, and derivative at sty.
#           The variables stp, f, g contain the values of the step,
#           function, and derivative at stp.
            stx = 0
            fx = finit
            gx = ginit
            sty = 0
            fy = finit
            gy = ginit
            stmin = 0
            stmax = stp + self.xtrapu * stp
            self.task = 'FG'
            self.save((stage, ginit, gtest, gx,
                       gy, finit, fx, fy, stx, sty,
                       stmin, stmax, width, width1))
            stp = self.determine_step(stp)
            #return stp, f, g
            return stp
        else:
            if self.isave[0] == 1:
                self.bracket = True
            else:
                self.bracket = False
            stage = self.isave[1]
            (ginit, gtest, gx, gy, finit, fx, fy, stx, sty, stmin, stmax, \
            width, width1) =self.dsave

#           If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
#           algorithm enters the second stage.
            ftest = finit + stp * gtest
            if stage == 1 and f < ftest and g >= 0.:
                stage = 2

#           Test for warnings.
            if self.bracket and (stp <= stmin or stp >= stmax):
                self.task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
            if self.bracket and stmax - stmin <= self.xtol * stmax:
                self.task = 'WARNING: XTOL TEST SATISFIED'
            if stp == self.stpmax and f <= ftest and g <= gtest:
                self.task = 'WARNING: STP = maxstep'
            if stp == self.stpmin and (f > ftest or g >= gtest):
                self.task = 'WARNING: STP = minstep'

#           Test for convergence.
            if f <= ftest and abs(g) <= c2 * (- ginit):
                self.task = 'CONVERGENCE'

#           Test for termination.
            if self.task[:4] == 'WARN' or self.task[:4] == 'CONV':
                self.save((stage, ginit, gtest, gx,
                           gy, finit, fx, fy, stx, sty,
                           stmin, stmax, width, width1))
                #return stp, f, g
                return stp

#              A modified function is used to predict the step during the
#              first stage if a lower function value has been obtained but
#              the decrease is not sufficient.
            #if stage == 1 and f <= fx and f > ftest:
#           #    Define the modified function and derivative values.
            #    fm =f - stp * gtest
            #    fxm = fx - stx * gtest
            #    fym = fy - sty * gtest
            #    gm = g - gtest
            #    gxm = gx - gtest
            #    gym = gy - gtest

#               Call step to update stx, sty, and to compute the new step.
            #    stx, sty, stp, gxm, fxm, gym, fym = self.update (stx, fxm, gxm, sty,
            #                                        fym, gym, stp, fm, gm,
            #                                        stmin, stmax)

#           #    Reset the function and derivative values for f.

            #    fx = fxm + stx * gtest
            #    fy = fym + sty * gtest
            #    gx = gxm + gtest
            #    gy = gym + gtest

            #else:
#           Call step to update stx, sty, and to compute the new step.

            stx, sty, stp, gx, fx, gy, fy= self.update(stx, fx, gx, sty,
                                               fy, gy, stp, f, g,
                                               stmin, stmax)


#           Decide if a bisection step is needed.

            if self.bracket:
                if abs(sty-stx) >= .66 * width1:
                    stp = stx + .5 * (sty - stx)
                width1 = width
                width = abs(sty - stx)

#           Set the minimum and maximum steps allowed for stp.

            if self.bracket:
                stmin = min(stx, sty)
                stmax = max(stx, sty)
            else:
                stmin = stp + self.xtrapl * (stp - stx)
                stmax = stp + self.xtrapu * (stp - stx)

#           Force the step to be within the bounds maxstep and minstep.

            stp = max(stp, self.stpmin)
            stp = min(stp, self.stpmax)

            if (stx == stp and stp == self.stpmax and stmin > self.stpmax):
                self.no_update = True
#           If further progress is not possible, let stp be the best
#           point obtained during the search.

            if (self.bracket and stp < stmin or stp >= stmax) \
               or (self.bracket and stmax - stmin < self.xtol * stmax):
                stp = stx

#           Obtain another function and derivative.

            self.task = 'FG'
            self.save((stage, ginit, gtest, gx,
                       gy, finit, fx, fy, stx, sty,
                       stmin, stmax, width, width1))
            return stp

    def update(self, stx, fx, gx, sty, fy, gy, stp, fp, gp,
               stpmin, stpmax):
        sign = gp * (gx / abs(gx))

#       First case: A higher function value. The minimum is bracketed.
#       If the cubic step is closer to stx than the quadratic step, the
#       cubic step is taken, otherwise the average of the cubic and
#       quadratic steps is taken.
        if fp > fx:  #case1
            self.case = 1
            theta = 3. * (fx - fp) / (stp - stx) + gx + gp
            s = max(abs(theta), abs(gx), abs(gp))
            gamma = s * np.sqrt((theta / s) ** 2. - (gx / s) * (gp / s))
            if stp < stx:
                gamma = -gamma
            p = (gamma - gx) + theta
            q = ((gamma - gx) + gamma) + gp
            r = p / q
            stpc = stx + r * (stp - stx)
            stpq = stx + ((gx / ((fx - fp) / (stp-stx) + gx)) / 2.) \
                   * (stp - stx)
            if (abs(stpc - stx) < abs(stpq - stx)):
               stpf = stpc
            else:
               stpf = stpc + (stpq - stpc) / 2.

            self.bracket = True

#       Second case: A lower function value and derivatives of opposite
#       sign. The minimum is bracketed. If the cubic step is farther from
#       stp than the secant step, the cubic step is taken, otherwise the
#       secant step is taken.

        elif sign < 0:  #case2
            self.case = 2
            theta = 3. * (fx - fp) / (stp - stx) + gx + gp
            s = max(abs(theta), abs(gx), abs(gp))
            gamma = s * np.sqrt((theta / s) ** 2 - (gx / s) * (gp / s))
            if stp > stx:
                 gamma = -gamma
            p = (gamma - gp) + theta
            q = ((gamma - gp) + gamma) + gx
            r = p / q
            stpc = stp + r * (stx - stp)
            stpq = stp + (gp / (gp - gx)) * (stx - stp)
            if (abs(stpc - stp) > abs(stpq - stp)):
               stpf = stpc
            else:
               stpf = stpq
            self.bracket = True

#       Third case: A lower function value, derivatives of the same sign,
#       and the magnitude of the derivative decreases.

        elif abs(gp) < abs(gx):  #case3
            self.case = 3
#           The cubic step is computed only if the cubic tends to infinity
#           in the direction of the step or if the minimum of the cubic
#           is beyond stp. Otherwise the cubic step is defined to be the
#           secant step.

            theta = 3. * (fx - fp) / (stp - stx) + gx + gp
            s = max(abs(theta), abs(gx), abs(gp))

#           The case gamma = 0 only arises if the cubic does not tend
#           to infinity in the direction of the step.

            gamma = s * np.sqrt(max(0.,(theta / s) ** 2-(gx / s) * (gp / s)))
            if stp > stx:
                gamma = -gamma
            p = (gamma - gp) + theta
            q = (gamma + (gx - gp)) + gamma
            r = p / q
            if r < 0. and gamma != 0:
               stpc = stp + r * (stx - stp)
            elif stp > stx:
               stpc = stpmax
            else:
               stpc = stpmin
            stpq = stp + (gp / (gp - gx)) * (stx - stp)

            if self.bracket:

#               A minimizer has been bracketed. If the cubic step is
#               closer to stp than the secant step, the cubic step is
#               taken, otherwise the secant step is taken.

                if abs(stpc - stp) < abs(stpq - stp):
                    stpf = stpc
                else:
                    stpf = stpq
                if stp > stx:
                    stpf = min(stp + .66 * (sty - stp), stpf)
                else:
                    stpf = max(stp + .66 * (sty - stp), stpf)
            else:

#               A minimizer has not been bracketed. If the cubic step is
#               farther from stp than the secant step, the cubic step is
#               taken, otherwise the secant step is taken.

                if abs(stpc - stp) > abs(stpq - stp):
                   stpf = stpc
                else:
                   stpf = stpq
                stpf = min(stpmax, stpf)
                stpf = max(stpmin, stpf)

#       Fourth case: A lower function value, derivatives of the same sign,
#       and the magnitude of the derivative does not decrease. If the
#       minimum is not bracketed, the step is either minstep or maxstep,
#       otherwise the cubic step is taken.

        else:  #case4
            self.case = 4
            if self.bracket:
                theta = 3. * (fp - fy) / (sty - stp) + gy + gp
                s = max(abs(theta), abs(gy), abs(gp))
                gamma = s * np.sqrt((theta / s) ** 2 - (gy / s) * (gp / s))
                if stp > sty:
                    gamma = -gamma
                p = (gamma - gp) + theta
                q = ((gamma - gp) + gamma) + gy
                r = p / q
                stpc = stp + r * (sty - stp)
                stpf = stpc
            elif stp > stx:
                stpf = stpmax
            else:
                stpf = stpmin

#       Update the interval which contains a minimizer.

        if fp > fx:
            sty = stp
            fy = fp
            gy = gp
        else:
            if sign < 0:
                sty = stx
                fy = fx
                gy = gx
            stx = stp
            fx = fp
            gx = gp
#       Compute the new step.

        stp = self.determine_step(stpf)

        return stx, sty, stp, gx, fx, gy, fy

    def determine_step(self, stp):
        dr = stp - self.old_stp
        if abs(pymax(self.pk) * dr) > self.maxstep:
            dr /= abs((pymax(self.pk) * dr) / self.maxstep)
        stp = self.old_stp + dr
        return stp

    def save(self, data):
        if self.bracket:
            self.isave[0] = 1
        else:
            self.isave[0] = 0
        self.isave[1] = data[0]
        self.dsave = data[1:]
