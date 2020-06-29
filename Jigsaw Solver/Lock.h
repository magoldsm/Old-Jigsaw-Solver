#pragma once

using Indices = std::vector<int>;

void Lock(const GTransform& g0, const Curve& CDelta, const Curve& CtildeDelta, GTransform& gLock, double K3, Indices& D_Delta_3_Ics, Indices& Dtilde_Delta_3_Ics, Indices& D_Delta_2_Ics, Indices& Dtilde_Delta_2_Ics, LRESULT& plotHandle);
