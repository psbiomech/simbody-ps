#ifndef SimTK_SIMBODY_SMOOTH_SPHERE_CLOSED_HALFSPACE_FORCE_H_
#define SimTK_SIMBODY_SMOOTH_SPHERE_CLOSED_HALFSPACE_FORCE_H_

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
#include "simbody/internal/Force.h"

namespace SimTK {

    class SmoothSphereClosedHalfSpaceForceImpl;

    /** This class models the forces generated by simple point contacts between a
    sphere and a half space. The contact model is a smooth (i.e., twice
    continuously differentiable) approximation of the HuntCrossleyForce already
    available in Simbody. The proposed implementation was designed for use with
    gradient-based optimization algorithms and algorithmic/automatic
    differentiation. To this aim, conditional if statements were approximated by
    using hyperbolic tangent functions. For example, the following if statement:
    <pre>     y = 0, if x < d </pre>
    <pre>     y = a, if x >= d </pre>
    can be approximated by:
    <pre>     f = 0.5 + 0.5 tanh(b(x-d)) </pre>
    <pre>     y = a f </pre>
    where b is a parameter that determines the smoothness of the transition.

    This force does not rely on a GeneralContactSubsystem. Instead, it assumes
    contact between a sphere and a half space, where the half space is defined as
    x > 0 (that is, all coordinates with positive x are in contact) in half space's
    frame. Commonly, the gravity acceleration vector is (0, -g, 0) (g > 0), and the
    half space represents ground. To achieve this, the half space should be rotated
    -90 degrees about the Z-Axis. This can be done by defining the half space frame
    as: @code Transform exampleHalfSpaceFrame(Rotation(-0.5*Pi, ZAxis), Vec3(0))
    @endcode

    The contact model includes components for the normal restoring force,
    dissipation in the material, and surface friction. The force is only applied to
    point contacts.

    To use it, do the following:

    <ol>
    <li>Add a GeneralForceSubsystem to a MultibodySystem.</li>
    <li>Add a SmoothSphereClosedHalfSpaceForce to the GeneralForceSubsystem.</li>
    <li>Add a MobilizedBody for the contact sphere and for the contact half
    space, and call setContactSphereBody(), setContactSphereLocationInBody(),
    setContactSphereRadius(), setContactHalfSpaceBody(),
    setContactHalfSpaceFrame(), and setParameters(). </li>
    </ol>

    <h1>Normal Force Components</h1>

    The normal restoring force (Hertz force) is given by:
    <pre>     fh_pos = (4/3) k (R k)^(1/2) ((x^2+cf)^(1/2))^(3/2) </pre>
    where k = 0.5 E^(2/3) with E the plain strain modulus, which is assumed
    identical for both contacting materials (i.e., sphere and half space), x is
    penetration depth, R is sphere radius, and cf (default is 1e-5) is a constant
    that enforces non-null derivatives. To smoothly transition between periods with
    and without sphere-half space contact, we use a tanh function:
    <pre>     fh_smooth = fh_pos (1/2+(1/2)tanh(bd x)) </pre>
    where bd (default is 300) is a parameter that determines the smoothness of the
    tanh transition. The graph below compares the smooth approximation with
    respect to the original model.

    \htmlonly <style>div.image img[src="SmoothSphereClosedHalfSpaceForce_HertzForce.png"]{width:750px;}</style> \endhtmlonly
    @image html SmoothSphereClosedHalfSpaceForce_HertzForce.png "Curves produced using E=1e6, R=0.8, cf=1e-5, and bd=300"

    The dissipation force is combined with the normal restoring force
    (Hunt-Crossley force) as follows:
    <pre>     fhc_pos = fh_smooth (1+(3/2) c v) </pre>
    where c is dissipation and v is penetration rate. To smoothly transition
    between null and positive Hunt-Crossley force, we use a tanh function:
    <pre>     fhc_smooth = fhc_pos (1/2+(1/2) tanh(bv (v+(2/(3 c))))) </pre>
    where bv (default is 50) is a parameter that determines the smoothness of the
    tanh transition. The graph below compares the smooth approximation with
    respect to the original model.

    \htmlonly <style>div.image img[src="SmoothSphereClosedHalfSpaceForce_HuntCrossleyForce.png"]{width:750px;}</style> \endhtmlonly
    @image html SmoothSphereClosedHalfSpaceForce_HuntCrossleyForce.png "Curves produced using x=0.1, E=1e6, R=0.8, c=2, cf=1e-5, bd=300, and bv=50"

    <h1>Friction Force</h1>

    The friction force is given by:
    <pre> ff = fhc_smooth [min(vs/vt,1) (ud+2(us-ud)/(1+(vs/vt)^2))+uv vs] </pre>
    where vs is the slip velocity of the two bodies at the contact point (see
    below), vt is a transition velocity (see below), and us, ud, and uv are the
    coefficients of static, dynamic, and viscous friction, respectively.

    The slip velocity is defined as the norm of the tangential velocity. To
    enforce non-null derivatives, we added the small positive constant cf (default
    is 1e-5):
    <pre> vs = sqrt(vtangent[1]^2 + vtangent[2]^2 + vtangent[3]^2 + cf) </pre>
    where vtangent is the tangential velocity.

    Because the friction force is a continuous function of the slip velocity, this
    model cannot represent stiction; as long as a tangential force is applied, the
    two bodies will move relative to each other. There will always be a nonzero
    drift, no matter how small the force is. The transition velocity vt acts as an
    upper limit on the drift velocity. By setting vt to a sufficiently small value,
    the drift velocity can be made arbitrarily small, at the cost of making the
    equations of motion very stiff.

    <h1>Contact Force</h1>

    The contact force is given by:
    <pre> force = fhc_smooth*normal + ff*vtangent/vs </pre>
    where normal determines the direction of the normal to the contact half space
    (the normal points in the direction of contact).

    <h1>Potential energy</h1>

    The potential energy is the integral of the normal restoring force (Hertz
    force). Due to the smooth approximation, there is no exact expression for the
    potential energy. We therefore made a comprise by providing an approximation
    for the potential energy. Specifically, we used the original expression but
    replaced the original Hertz force by our smooth approximation as follows:
    <pre> pe = (2/5)*fh_smooth*x </pre>
    This approximation results in the error depicted in the graph below when
    differentiating the potential energy with respect to the penetration, the user
    should thus be aware of the limitations of this model that was primarily
    designed for optimization problems.

    \htmlonly <style>div.image img[src="SmoothSphereClosedHalfSpaceForce_HertzForceEnergyError.png"]{width:750px;}</style> \endhtmlonly
    @image html SmoothSphereClosedHalfSpaceForce_HertzForceEnergyError.png "Curves produced using E=1e6, R=0.8, cf=1e-5, and bd=300"
    */
    class SimTK_SIMBODY_EXPORT SmoothSphereClosedHalfSpaceForce : public Force {
    public:
        /** Create a smooth sphere to half space Hunt-Crossley contact model.
        @param forces the subsystem that will own this
        SmoothSphereClosedHalfSpaceForce element */
        SmoothSphereClosedHalfSpaceForce(GeneralForceSubsystem& forces);
        /** Set the contact material parameters.
        @param stiffness the stiffness constant (i.e., plain strain modulus),
            default is 1 (N/m^2)
        @param dissipation the dissipation coefficient, default is 0 (s/m)
        @param staticFriction the coefficient of static friction, default is 0
        @param dynamicFriction the coefficient of dynamic friction, default is 0
        @param viscousFriction the coefficient of viscous friction, default is 0
        @param transitionVelocity the transition velocity, default is 0.01 (m/s)
        @param cf the constant that enforces non-null derivatives, default is 1e-5
        @param bd the parameter that determines the smoothness of the transition of
            the tanh used to smooth the Hertz force, default is 300
        @param bv the parameter that determines the smoothness of the transition of
            the tanh used to smooth the Hunt-Crossley force, default is 50 */
        void setParameters(Real stiffness, Real dissipation, Real staticFriction,
            Real dynamicFriction, Real viscousFriction, Real transitionVelocity,
            Real cf, Real bd, Real bv, Real fpTanhCoeff,
            Vec3 fpCenter, Real fpXdim, Real fpZdim);
        /** Set the stiffness constant (i.e., plain strain modulus), default is 1
            (N/m^2). */
        void setStiffness(Real stiffness);
        /** Set the dissipation coefficient, default is 0 (s/m). */
        void setDissipation(Real dissipation);
        /** Set the coefficient of static friction, default is 0. */
        void setStaticFriction(Real staticFriction);
        /** Set the coefficient of dynamic friction, default is 0. */
        void setDynamicFriction(Real dynamicFriction);
        /** Set the coefficient of viscous friction, default is 0. */
        void setViscousFriction(Real viscousFriction);
        /** Set the transition velocity, default is 0.01 (m/s). */
        void setTransitionVelocity(Real transitionVelocity);
        /** Set the constant that enforces non-null derivatives, default is 1e-5.*/
        void setConstantContactForce(Real cf);
        /** Set the parameter that determines the smoothness of the transition
            of the tanh used to smooth the Hertz force. The larger the steeper the
            transition but also the more discontinuous, default is 300. */
        void setHertzSmoothing(Real bd);
        /** Set the parameter that determines the smoothness of the transition
            of the tanh used to smooth the Hunt-Crossley force. The larger the
            steeper the transition but also the more discontinuous, default
            is 50. */
        void setHuntCrossleySmoothing(Real bv);
        // Set the coefficient of the tanh boundary function arguments
        void setFPTanhCoeff(Real fpTanhCoeff);
        /** Set the force platform location and dimensions **/
        void setFPInfo(Vec3 fpCenter, Real fpXdim, Real fpZdim);
        /** Set the MobilizedBody to which the contact sphere is attached. */
        void setContactSphereBody(MobilizedBody bodyInput1);
        /** Set the location of the contact sphere in the body frame. */
        void setContactSphereLocationInBody(Vec3 locationSphere);
        /** Set the radius of the contact sphere. */
        void setContactSphereRadius(Real radius);
        /** Set the MobilizedBody to which the contact half space is attached. */
        void setContactHalfSpaceBody(MobilizedBody bodyInput2);
        /** Set the transform of the contact half space in the body frame. */
        void setContactHalfSpaceFrame(Transform halfSpaceFrame);
        /** Get the MobilizedBody to which the contact sphere is attached. */
        MobilizedBody getBodySphere();
        /** Get the MobilizedBody to which the contact half space is attached. */
        MobilizedBody getBodyHalfSpace();
        /** Get the location of the contact sphere in the body frame. */
        Vec3 getContactSphereLocationInBody();
        /** Get the radius of the contact sphere. */
        Real getContactSphereRadius();
        /** Get the transform of the contact half space. */
        Transform getContactHalfSpaceTransform();
        SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(SmoothSphereClosedHalfSpaceForce,
            SmoothSphereClosedHalfSpaceForceImpl, Force);
    };

} // namespace SimTK

#endif // SimTK_SIMBODY_SMOOTH_SPHERE_CLOSED_HALFSPACE_FORCE_H_

