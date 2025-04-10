#include "../HyRRT.h"
#include "ompl/base/StateSpace.h"
#include "../src/HyRRT.cpp"
#include "ompl/control/spaces/RealVectorControlSpace.h"
#include <fstream>
#include "../src/HybridStateSpace.cpp"
#include <ompl/control/ODESolver.h>

/** \brief Function computes the Pythagorean distance between two states. */
double distanceFunc(ompl::base::State *state1, ompl::base::State *state2)
{
    double dist = 0;
    dist = sqrt(pow(state1->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0] - state2->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0], 2) + pow(state1->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1] - state2->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1], 2));
    return fabs(dist);
}

/** \brief Jump set is true whenever the ball is on or below the surface and has a downwards velocity. */
bool jumpSet(ompl::control::HyRRT::Motion *motion)
{
    auto *motion_state = motion->state->as<ompl::base::CompoundState>()->as<ompl::base::RealVectorStateSpace::StateType>(0);
    double velocity = motion_state->values[1];
    double pos_cur = motion_state->values[0];

    if (pos_cur <= 0 && velocity <= 0)
        return true;
    else
        return false;
}

/** \brief Flow set is true whenever the ball is above the surface or has an upwards velocity. */
bool flowSet(ompl::control::HyRRT::Motion *motion)
{
    return !jumpSet(motion);
}

/** \brief Unsafe set is true whenever the ball is above 10 units from the ground, to reduce time spent planning. */
bool unsafeSet(ompl::control::HyRRT::Motion *motion)
{
    double u = motion->control->as<ompl::control::RealVectorControlSpace::ControlType>()->values[0];
    if (u > 5 || u < 0)
        return true;
    return false;
}

/** \brief Simulates the dynamics of the ball when in flow regime, with no nonnegligble forces other than gravity. */
ompl::base::State *continuousSimulator(std::vector<double> inputs, ompl::base::State *x_cur, double tFlow, ompl::base::State *new_state)
{
    double velocity = x_cur->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1];
    double acceleration_cur = x_cur->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2];
    double pos_cur = x_cur->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0];

    double v = velocity + (acceleration_cur) * tFlow;                               // v = v0 + at
    double x = pos_cur + velocity * tFlow + (acceleration_cur) * pow(tFlow, 2) / 2; // x = v0 * t + 1/2(at^2)

    new_state->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0] = x;
    new_state->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1] = v;
    new_state->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2] = acceleration_cur;

    return new_state;
}

void flowODE(const ompl::control::ODESolver::StateType& q, const ompl::control::Control* c, ompl::control::ODESolver::StateType& qdot)
{
    // Ensure qdot is the same size as q.  Zero out all values.
    qdot.resize(q.size(), 0);
 
    qdot[0] = q[1];            // x-dot
    qdot[1] = q[2];            // y-dot
    qdot[2] = 0;    // theta-dot
}

/** \brief Simulates the dynamics of the ball when in jump regime, with input from the surface. */
ompl::base::State *discreteSimulator(ompl::base::State *x_cur, const ompl::control::Control *u, ompl::base::State *new_state)
{
    double velocity = -0.8 * x_cur->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1];

    new_state->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0] = x_cur->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0];
    new_state->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1] = velocity - u->as<ompl::control::RealVectorControlSpace::ControlType>()->values[0];
    new_state->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2] = x_cur->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2];
    return new_state;
}


int main()
{
    std::uint_fast32_t seed = 1;
    ompl::RNG::setSeed(seed);
    // Set the bounds of space
    ompl::base::RealVectorStateSpace *statespace = new ompl::base::RealVectorStateSpace(0);
    statespace->addDimension(-10, 10);
    statespace->addDimension(-100, 100);
    statespace->addDimension(-10, 5);
    ompl::base::StateSpacePtr stateSpacePtr(statespace);

    ompl::base::HybridStateSpace *hybridSpace = new ompl::base::HybridStateSpace(stateSpacePtr);
    ompl::base::StateSpacePtr hybridSpacePtr(hybridSpace);

    // Define control space
    ompl::control::RealVectorControlSpace *flowControlSpace = new ompl::control::RealVectorControlSpace(hybridSpacePtr, 1);
    ompl::control::RealVectorControlSpace *jumpControlSpace = new ompl::control::RealVectorControlSpace(hybridSpacePtr, 1);

    ompl::base::RealVectorBounds flowBounds(1);
    flowBounds.setLow(0, 0);
    flowBounds.setHigh(0, 0);
    flowControlSpace->setBounds(flowBounds);

    ompl::base::RealVectorBounds jumpBounds(1);
    jumpBounds.setLow(0, 0);
    jumpBounds.setHigh(0, 5);
    jumpControlSpace->setBounds(jumpBounds);

    ompl::control::RealVectorControlUniformSampler flowControlSampler(flowControlSpace);
    flowControlSpace->setControlSamplerAllocator([flowControlSpace](const ompl::control::ControlSpace *space) -> ompl::control::ControlSamplerPtr {
        return std::make_shared<ompl::control::RealVectorControlUniformSampler>(space);
    });

    ompl::control::RealVectorControlUniformSampler jumpControlSampler(jumpControlSpace);     // Doesn't do anything because the bounds for jump input are just [0, 0], but here for demonstration
    jumpControlSpace->setControlSamplerAllocator([jumpControlSpace](const ompl::control::ControlSpace *space) -> ompl::control::ControlSamplerPtr {
        return std::make_shared<ompl::control::RealVectorControlUniformSampler>(space);
    });

    ompl::control::ControlSpacePtr flowControlSpacePtr(flowControlSpace);
    ompl::control::ControlSpacePtr jumpControlSpacePtr(jumpControlSpace);

    ompl::control::CompoundControlSpace *controlSpace = new ompl::control::CompoundControlSpace(hybridSpacePtr);
    controlSpace->addSubspace(flowControlSpacePtr);
    controlSpace->addSubspace(jumpControlSpacePtr);

    ompl::control::ControlSpacePtr controlSpacePtr(controlSpace);

    // Construct a space information instance for this state space
    ompl::control::SpaceInformationPtr si(new ompl::control::SpaceInformation(hybridSpacePtr, controlSpacePtr));
    ompl::control::ODESolverPtr odeSolver (new ompl::control::ODEBasicSolver<> (si, &flowODE));
    
    si->setStatePropagator(ompl::control::ODESolver::getStatePropagator(odeSolver));
    si->setPropagationStepSize(0.01);
    si->setup();
    

    // // Construct a space information instance for this state space
    // ompl::base::SpaceInformationPtr si(new ompl::base::SpaceInformation(hybridSpacePtr));

    // si->setup();

    // Set start state to be 1 unit above the floor with zero velocity and gravitaional acceleration of 9.81
    ompl::base::ScopedState<> start(hybridSpacePtr);
    start->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0] = 1;
    start->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1] = 0;
    start->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[2] = -9.81;

    // Set goal state to be on floor with a zero velocity
    ompl::base::ScopedState<> goal(hybridSpacePtr);
    goal->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[0] = 0;
    goal->as<ompl::base::HybridStateSpace::StateType>()->as<ompl::base::RealVectorStateSpace::StateType>(0)->values[1] = 0;

    // Create a problem instance
    ompl::base::ProblemDefinitionPtr pdef(new ompl::base::ProblemDefinition(si));

    // Set the start and goal states
    pdef->setStartAndGoalStates(start, goal);

    ompl::control::HyRRT cHyRRT(si);

    // Set the problem instance for our planner to solve
    cHyRRT.setProblemDefinition(pdef);
    cHyRRT.setup();

    // Set parameters
    cHyRRT.setDiscreteSimulator(discreteSimulator);
    cHyRRT.setDistanceFunction(distanceFunc);
    cHyRRT.setFlowSet(flowSet);
    cHyRRT.setJumpSet(jumpSet);
    cHyRRT.setTm(0.5);
    cHyRRT.setFlowStepDuration(0.01);
    cHyRRT.setFlowInputRange(std::vector<double>{0}, std::vector<double>{0});   // If input is a single value, only that value will every be used
    cHyRRT.setJumpInputRange(std::vector<double>{0}, std::vector<double>{5});
    cHyRRT.setUnsafeSet(unsafeSet);

    // attempt to solve the planning problem within 10 seconds
    ompl::base::PlannerStatus solved = cHyRRT.solve(ompl::base::timedPlannerTerminationCondition(100));
    // print path to RViz2 data file
    std::ofstream outFile("../../examples/visualize/src/points.txt");
    pdef->getSolutionPath()->as<ompl::control::PathControl>()->printAsMatrix(outFile);
    cout << "solution status: " << solved << endl;
}