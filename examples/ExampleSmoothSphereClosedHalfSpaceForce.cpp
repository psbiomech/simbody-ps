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

#include "Simbody.h"

#include <cmath>

using namespace SimTK;
using namespace std;


// Report the contact forces
class MyForceReporter : public PeriodicEventReporter {
public:
    MyForceReporter(const MultibodySystem& system, Real period)
        : PeriodicEventReporter(period), system(system) {}

    virtual void handleEvent(const State& state) const override {

        cout << "\nt=" << state.getTime() << endl;
        for (int b = 1; b <= 9; ++b) {
            cout << "\tF[" << b << "]=" << Vec3(system.getRigidBodyForces(state, Stage::Dynamics)[b][1]) << endl;
        }
    }
private:
    const MultibodySystem& system;
};


int main() {

    // Create the universe
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    const Vec3 gravity = Vec3(0, -9.8, 0);
    Force::UniformGravity(forces, matter, gravity, 0);

    // Define properties of the contact sphere and half spaces
    const Real radius = 0.1;
    const Real k = 1E6;
    const Real stiffness = 0.5*std::pow(k, 2.0/3.0);
    const Real dissipation = 0.5;
    const Real us = 1.0;
    const Real ud = 0.5;
    const Real uv = 0.1;
    const Real vt = 0.001;
    const Real cf = 1e-5;
    const Real bd = 300;
    const Real bv = 50;

    // Define the diagonal corners of the closed half space and the
    // coefficient of the tanh argument
    const Vec3 fpC1(-0.5, 0.5, -0.5);
    const Vec3 fpC2(0.5, 1, 0.5);
    const Real fpcoeff = 1E6;

    // Platform dimensions, orientation and centre
    Vec3 fpcentre = fpC1 + (fpC2 - fpC1) / 2;
    Real fpxlen = sqrt(pow(abs(fpC2[0] - fpC1[0]), 2) + pow(abs(fpC2[1] - fpC1[1]), 2));
    Real fpzlen = sqrt(pow(abs(fpC2[2] - fpC1[2]), 2) + pow(abs(fpC2[1] - fpC1[1]), 2));
    Vec3 fpdims(fpxlen, 0, fpzlen);
    Real fpxrot = atan((fpC2[1] - fpC1[1]) / (fpC2[2] - fpC1[2]));
    Real fpzrot = atan((fpC2[1] - fpC1[1]) / (fpC2[0] - fpC1[0]));
    Vec3 fpangles(fpxrot, 0, fpzrot);
    Rotation fporient;
    fporient.setRotationFromAngleAboutX(fpxrot);
    fporient.setRotationFromAngleAboutY(0);
    fporient.setRotationFromAngleAboutZ(fpzrot);


    // Create several bodies and attach contact spheres
    const int nspheres = 9;
    Body::Rigid bodies[nspheres];
    MobilizedBody::Translation spheres[nspheres];
    for (int b = 0; b < nspheres; ++b) {
        Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
        MobilizedBody::Translation sphere(matter.updGround(),
            Transform(), body, Transform());
        sphere.addBodyDecoration(Transform(), DecorativeSphere(radius));
        bodies[b] = body;
        spheres[b] = sphere;
    }
    

    // Create a body and attach the closed half space
    Body::Rigid bodyhs1(MassProperties(0, Vec3(0, 0, 0), Inertia(0)));
    MobilizedBody::Weld halfspace1(matter.updGround(),
        Transform(fporient, fpcentre), bodyhs1, Transform());
    halfspace1.addBodyDecoration(Transform(),
        DecorativeBrick(fpdims / 2).setColor(Red));

    // Create a body and attach the open half space
    Body::Rigid bodyhs2(MassProperties(0, Vec3(0), Inertia(0)));
    MobilizedBody::Weld halfspace2(matter.updGround(),
        Transform(), bodyhs2, Transform());
    //halfspace2.addBodyDecoration(Transform(),
    //    DecorativeBrick(Vec3(5, 0, 5)).setColor(Yellow));


    // Model contact between the spheres and closed half space using smooth 
    // Hunt-Crossley forces
    for (int b = 0; b < nspheres; ++b) {
        SmoothSphereClosedHalfSpaceForce hc_smooth(forces);
        hc_smooth.setParameters(k, dissipation, us, ud, uv, vt, cf, bd, bv,
            fpcoeff, fpC1, fpC2);
        hc_smooth.setContactSphereBody(spheres[b]);
        hc_smooth.setContactSphereLocationInBody(Vec3(0));
        hc_smooth.setContactSphereRadius(radius);
        hc_smooth.setContactHalfSpaceFrame(Transform(
            Rotation(-0.5 * Pi, ZAxis), Vec3(0)));
        hc_smooth.setContactHalfSpaceBody(halfspace1);
    }
 
    // Model contact between the spheres and open half space using smooth 
    // Hunt-Crossley forces
    for (int b = 0; b < nspheres; ++b) {
        SmoothSphereHalfSpaceForce hc_smooth(forces);
        hc_smooth.setParameters(k, dissipation, us, ud, uv, vt, cf, bd, bv);
        hc_smooth.setContactSphereBody(spheres[b]);
        hc_smooth.setContactSphereLocationInBody(Vec3(0));
        hc_smooth.setContactSphereRadius(radius);
        hc_smooth.setContactHalfSpaceFrame(Transform(
            Rotation(-0.5 * Pi, ZAxis), Vec3(0)));
        hc_smooth.setContactHalfSpaceBody(halfspace2);
    }

    // Create a visualiser, also report contact forces
    Visualizer viz(system);
    viz.setBackgroundType(Visualizer::BackgroundType(Visualizer::GroundAndSky));
    //viz.setBackgroundColor(White);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));
    system.addEventReporter(new MyForceReporter(system, 0.01));

    // Construct the system and return the default state
    State state = system.realizeTopology();
    system.realizeModel(state);

    // Initialise the state of the sphere
    std::array<double, nspheres> init_x
        = { -1., -0.5, 0, 0, 0, 0, 0, 0.5, 1. };
    std::array<double, nspheres> init_z 
        = { 0, 0, -0.5, -1., 0, 0.5, 1., 0, 0 };
    double init_y = 2;
    for (int b = 0; b < nspheres; ++b) {
        Vec3 init_pos = { init_x[b], init_y, init_z[b] };
        spheres[b].setQFromVector(state, Vector(init_pos));
    }
    
    // Simulate and report forces
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(Real(1e-3));
    TimeStepper ts(system, integ);
    ts.initialize(state); // set IC's
    ts.stepTo(5.0);

}

