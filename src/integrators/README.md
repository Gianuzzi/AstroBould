# Integrators Module

This module provides a unified interface to a large collection of numerical integrators for ODE systems, including explicit Runge–Kutta methods, implicit schemes, symplectic leapfrog, Bulirsch–Stoer, and adaptive embedded RK pairs.

The integrator is selected by an integer parameter **`i`**, and the module automatically dispatches to the correct algorithm depending on the value of `i`.

The design goal is to allow a single simulation code to choose integrators dynamically without changing the calling sequence.

---

## Integrator Selection

The rule for selecting the integrator is:

* **`i = 0`** → full Bulirsch–Stoer (BS)
* **`i = –1, –2`** → Leapfrog variants
* **`i = –3`** → BS with Simplified Midpoint (fixed sequence)
* **`i ≤ –4`** → explicit/implicit Runge–Kutta methods (fixed step)
* **`i ≥ 1`** → embedded Runge–Kutta pairs (adaptive step)

Formally:

```fortran
select case (i)

case (0)
    ! Bulirsch–Stoer (adaptive, extrapolated)

case (-1)
    ! Leapfrog drift–kick–drift (2nd order)

case (-2)
    ! Leapfrog kick–drift–kick (2nd order)

case (-3)
    ! BS simplified-midpoint (extrapolation range)

case (: -4)
    ! Explicit/implicit RK method number abs(i+3)

case (1 :)
    ! Embedded RK pair number i
end select
```

Both Leapfrog and RK methods can optionally be executed with an adaptive timestep using the provided error-control driver.
Adaptive stepping **does not change the formal order**, but removes symplecticity for the Leapfrog schemes.

---

## List of Available Integrators

The following table lists all integrators exposed by the module, with their formal order.
Symplecticity applies only for constant timestep.

| Option | Integrator                         | Order                           |
| ------ | ---------------------------------- | ------------------------------- |
| -29    | Abbas6                             | 6                               |
| -28    | RK6                                | 6                               |
| -27    | Gauss–Legendre6 (Implicit)         | 6                               |
| -26    | Nyström5                           | 5                               |
| -25    | RK_four_oct4                       | 4                               |
| -24    | Gauss–Legendre4 (Implicit)         | 4                               |
| -23    | RK4                                | 4                               |
| -22    | Lobatto4 (Implicit)                | 4                               |
| -21    | Ralston4 (Implicit)                | 4                               |
| -20    | RK_implicit3 (Implicit)            | 3                               |
| -19    | Crouzeix3 (Implicit)               | 3                               |
| -18    | SSPRrk3                            | 3                               |
| -17    | Ralston3                           | 3                               |
| -16    | Heun3                              | 3                               |
| -15    | Runge_Kutta3                       | 3                               |
| -14    | Qin_Zhang2                         | 2                               |
| -13    | Kraaijevanger_Spijker2             | 2                               |
| -12    | Hammer_Hollingsworth2              | 2                               |
| -11    | Ralston2                           | 2                               |
| -10    | strange2                           | 2                               |
| -9     | midpoint2                          | 2                               |
| -8     | Heun2                              | 2                               |
| -7     | Crank–Nicolson2 (Implicit)         | 2                               |
| -6     | Euler_center2 (Implicit)           | 2                               |
| -5     | Euler_back1 (Implicit)             | 1                               |
| -4     | Euler1                             | 1                               |
| -3     | Bulirsch–Stoer Simplified Midpoint | 4–16 (extrapolated)             |
| -2     | Leapfrog (kick–drift–kick)         | 2 (symplectic only if dt fixed) |
| -1     | Leapfrog (drift–kick–drift)        | 2 (symplectic only if dt fixed) |
| 0      | Bulirsch–Stoer (full)              | 4–16 (extrapolated)             |
| 1      | Fehlberg1_2                        | (1, 2)                          |
| 2      | Heun_Euler2_1                      | (2, 1)                          |
| 3      | Fehlberg2_1                        | (2, 1)                          |
| 4      | Bogacki–Shampine3_2                | (3, 2)                          |
| 5      | Zonneveld4_3                       | (4, 3)                          |
| 6      | Merson4_3                          | (4, 3)                          |
| 7      | Fehlberg4_5                        | (4, 5)                          |
| 8      | Cash–Karp5_4                       | (5, 4)                          |
| 9      | Dormand–Prince5_4                  | (5, 4)                          |
| 10     | Verner6_5                          | (6, 5)                          |
| 11     | Fehlberg7_8                        | (7, 8)                          |
| 12     | Dormand–Prince8_7                  | (8, 7)                          |

---

## Adaptive Timestepping

All integrators may be used with the adaptive driver provided in the module.
Important notes:

### Leapfrog

* Formal order stays **2**
* **Not symplectic** when `dt` varies
* Energy conservation becomes similar to 2nd-order RK

### Embedded RK pairs

* Designed for adaptive stepping
* Order = first number in the pair
* Low-order solution provides error estimate

### Bulirsch–Stoer

* Always adaptive
* Extrapolated order increases with the number of midpoint steps
* Effective order typically between 4 and 16

---

## Practical Recommendations

* For **long-term Hamiltonian integrations**:
  use **Leapfrog with fixed dt**, or an explicit RK for non-symplectic scenarios.

* For **high accuracy, smooth problems**:
  use **Bulirsch–Stoer (0)** or **Verner6_5 (10)**.

* For **fast, robust, general-purpose ODE solving**:
  use **Dormand–Prince5_4 (9)**.

* For **stiff or semi-stiff problems**:
  use an **implicit Gauss–Legendre** or **Crank–Nicolson** scheme.

---

## Interface Overview

The module exposes one unified driver:

```
call integrate(t, y, dt_adap, dydt, dt, ynew, check_fun)
```
---