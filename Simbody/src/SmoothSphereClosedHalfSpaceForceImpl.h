#ifndef SimTK_SIMBODY_SMOOTH_SPHERE_CLOSED_HALFSPACE_FORCE_IMPL_H_
#define SimTK_SIMBODY_SMOOTH_SPHERE_CLOSED_HALFSPACE_FORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-19 Stanford University and the Authors.        *
 * Authors: Prasanna Sritharan, Antoine Falisse, Gil Serrancoli               *
 * Contributors: Peter Eastman                                                *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SmoothSphereClosedHalfSpaceForce.h"
#include "ForceImpl.h"

namespace SimTK {

    class SmoothSphereClosedHalfSpaceForceImpl : public ForceImpl {
    public:
        class Parameters {
        public:
            Parameters() : stiffness(1), dissipation(0), staticFriction(0),
                dynamicFriction(0), viscousFriction(0), transitionVelocity(0.01),
                cf(1e-5), bd(300), bv(50), fpTanhCoeff(1E6), 
                fpCorner1(Vec3(-0.5, 0, 0.5)), fpCorner2(Vec3(0.5, 0, -0.5)) {
            }
            Parameters(Real stiffness, Real dissipation, Real staticFriction,
                Real dynamicFriction, Real viscousFriction,
                Real transitionVelocity, Real cf, Real bd, Real bv, Real fpTanhCoeff,
                Vec3 fpCorner1, Vec3 fpCorner2) :
                stiffness(stiffness), dissipation(dissipation),
                staticFriction(staticFriction), dynamicFriction(dynamicFriction),
                viscousFriction(viscousFriction),
                transitionVelocity(transitionVelocity), cf(cf), bd(bd), bv(bv),
                fpTanhCoeff(fpTanhCoeff), fpCorner1(fpCorner1), 
                fpCorner2(fpCorner2) {
            }
            Real stiffness, dissipation, staticFriction, dynamicFriction,
                viscousFriction, transitionVelocity, cf, bd, bv, fpTanhCoeff;
            Vec3 fpCorner1, fpCorner2;
        };

        Real            contactSphereRadius;
        Vec3            contactSphereLocation;
        Transform       contactHalfSpaceFrame;
        MobilizedBody   bodySphere;
        MobilizedBody   bodyHalfSpace;
        Parameters      parameters;


        SmoothSphereClosedHalfSpaceForceImpl(GeneralForceSubsystem& subsystem);

        SmoothSphereClosedHalfSpaceForceImpl* clone() const override {
            return new SmoothSphereClosedHalfSpaceForceImpl(*this);
        }
        // Set the contact material parameters.
        void setParameters(Real stiffness, Real dissipation, Real staticFriction,
            Real dynamicFriction, Real viscousFriction, Real transitionVelocity,
            Real cf, Real bd, Real bv, Real fpTanhCoeff, 
            Vec3 fpCorner1, Vec3 fpCorner2);
        // Get parameters.
        const Parameters& getParameters() const;
        // Update parameters.
        Parameters& updParameters();
        // Set the stiffness constant (i.e., plain strain modulus), default is 1
        // N/m^2.
        void setStiffness(Real stiffness);
        // Set the dissipation coefficient, default is 0 s/m.
        void setDissipation(Real dissipation);
        // Set the coefficient of static friction, default is 0.
        void setStaticFriction(Real staticFriction);
        // Set the coefficient of dynamic friction, default is 0.
        void setDynamicFriction(Real dynamicFriction);
        // Set the coefficient of viscous friction, default is 0.
        void setViscousFriction(Real viscousFriction);
        // Set the transition velocity, default is 0.01 m/s.
        void setTransitionVelocity(Real transitionVelocity);
        // Set the constant that enforces non-null derivatives, default is 1e-5.
        void setConstantContactForce(Real cf);
        // Set the parameter that determines the smoothness of the transition
        // of the tanh used to smooth the Hertz force. The larger the steeper the
        // transition but also the more discontinuous, default is 300.
        void setHertzSmoothing(Real bd);
        // Set the parameter that determines the smoothness of the transition of
        // the tanh used to smooth the Hunt-Crossley force. The larger the steeper
        // the transition but also the more discontinuous, default is 50.
        void setHuntCrossleySmoothing(Real bv);
        // Set the coefficient of the tanh boundary function arguments
        void setFPTanhCoeff(Real fpTanhCoeff);
        // Set the diagnoal corners of the platform
        void setFPDiagonalCorners(Vec3 fpCorner1, Vec3 fpCorner2);
        // Set the MobilizedBody to which the contact sphere is attached.
        void setContactSphereBody(MobilizedBody bodyInput1);
        // Set the location of the contact sphere in the body frame.
        void setContactSphereLocationInBody(Vec3 locationContactSphere);
        // Set the radius of the contact sphere.
        void setContactSphereRadius(Real radius);
        // Set the MobilizedBody to which the contact half space is attached.
        void setContactHalfSpaceBody(MobilizedBody bodyInput2);
        // Set the transform of the contact half space in the body frame.
        void setContactHalfSpaceFrame(Transform halfSpaceFrame);
        // Get the MobilizedBody to which the contact sphere is attached.
        MobilizedBody getBodySphere();
        // Get the MobilizedBody to which the contact half space is attached.
        MobilizedBody getBodyHalfSpace();
        // Get the location of the contact sphere in the body frame.
        Vec3 getContactSphereLocationInBody();
        // Get the radius of the contact sphere.
        Real getContactSphereRadius();
        // Get the transform of the contact half space.
        Transform getContactHalfSpaceTransform();
        // Get the normal to the contact half space.
        void getNormalContactHalfSpace(const State& state,
            UnitVec3& normalContactHalfSpace) const;
        // Get the location of the contact sphere origin in the ground frame.
        void getContactSphereOrigin(const State& state, Vec3& contactPointPos)const;
        // Calculate contact force.
        void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const override;
        // Calculate potential energy.
        Real calcPotentialEnergy(const State& state) const override;
        void realizeTopology(State& state) const override;
    private:
        const GeneralForceSubsystem& subsystem;
        mutable CacheEntryIndex               energyCacheIndex;
    };

} // namespace SimTK

#endif // SimTK_SIMBODY_SMOOTH_SPHERE_CLOSED_HALFSPACE_FORCE_IMPL_H_

