/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2025, University of Santa Cruz Hybrid Systems Laboratory
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Author: Beverly Xu */

#ifndef OMPL_CONTROL_PLANNERS_RRT_HYRRT_
#define OMPL_CONTROL_PLANNERS_RRT_HYRRT_

#include "ompl/base/spaces/RealVectorStateSpace.h"
#include "ompl/control/ControlSpace.h"
#include "ompl/control/spaces/RealVectorControlSpace.h"
#include "ompl/datastructures/NearestNeighbors.h"
#include "ompl/control/planners/PlannerIncludes.h"
#include "ompl/control/Control.h"
#include "HybridStateSpace.h"

using namespace std;

namespace ompl
{
    namespace control
    {
        /**
           @anchor HyRRT
           @par Hybrid RRT (HyRRT) is an RRT algorithm that solves separate optimization
           problems associated with the flows and jumps of the systems, to solve a
           variety of robotic motion planning problems. As an RRT algorithm, HyRRT is
           probabilistically-complete. The logical flow of the algorithm is as follows:
           1. Initialize the tree with a start state.
           2. While a solution has not been found:
                a. Sample a random state.
                b. Find the nearest state in the tree to the random state.
                c. Extend from that state under either the flow or jump regimes, using
           the continuous or discrete simulators, respectively. d. Continue until the
           state is in collision, or the maximum flow time has been exceeded. e. If the
           state is not within the unsafe set, add the state to the tree. If the state
                   is in collision, proceed to jump.
            3. Return the solution.

           @par External documentation: https://ieeexplore.ieee.org/document/9992444
           DOI: [10.1109/CDC51059.2022.9992444]
        */

        /** \brief Hybrid Rapidly-exploring Random Trees */
        class HyRRT : public base::Planner
        {
        public:
            /** \brief Constructor */
            HyRRT(const control::SpaceInformationPtr &si);

            /** \brief Destructor */
            ~HyRRT() override;

            /** \brief Clear all allocated memory. */
            void clear() override;

            /** \brief Set the problem instance to solve */
            void setup() override;

            /** 
             * \brief Get the PlannerData object associated with this planner
             * @param data the PlannerData object storing the edges and vertices of the solution
             */
            void getPlannerData(base::PlannerData &data) const override;
            base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc) override;

            /**
             * \brief Set of available random input methods
             * @par See https://ompl.kavrakilab.org/classompl_1_1RNG.html for more information about each distribution
             */
            enum inputSamplingMethods_
            {
                UNIFORM_01,
                UNIFORM_REAL,
                UNIFORM_INT,
                GAUSSIAN_01,
                GAUSSIAN_REAL,
                HALF_NORMAL_REAL,
                HALF_NORMAL_INT,
                QUATERNION,
                EURLER_RPY
            };

            /// \brief Representation of a motion in the search tree
            class Motion
            {
            public:
                /// \brief Default constructor
                Motion() = default;

                /// \brief Constructor that allocates memory for the state
                Motion(const control::SpaceInformation *si) : state(si->allocState()), control(si->allocControl()) {}

                /// \brief Destructor
                ~Motion() = default;

                /// \brief The state contained by the motion
                base::State *state{nullptr};

                /// \brief The parent motion in the exploration tree
                Motion *parent{nullptr};

                /// \brief Pointer to the root of the tree this motion is
                /// contained in.
                const base::State *root{nullptr};

                /// \brief The integration steps defining the solution pair of the motion, between the parent and child vertices
                std::vector<base::State *> *solutionPair{nullptr};

                /** \brief The control contained by the motion */
                control::Control *control{nullptr};

                // /// \brief The inputs associated with the solution pair
                // std::vector<ompl::control::Control *> *inputs = new std::vector<ompl::control::Control *>();
            };

            /** \brief Free the memory allocated by this planner */
            void freeMemory();

            /** 
             * \brief Define the continuous interval over which the flow input is sampled. 
             * @param min the minimum input values
             * @param max the maximum input values
             */
            void setFlowInputRange(std::vector<double> min, std::vector<double> max)
            {
                int size = max.size() > min.size() ? max.size() : min.size();
                if (min.size() != max.size())
                    throw Exception("Max input value (maxFlowInputValue) and min input value (minFlowInputValue) must be of the same size");
                for (int i = 0; i < size; i++)
                {
                    if (min.at(0) > max.at(0))
                        throw Exception("Max input value must be greater than or equal to min input value");
                }
                minFlowInputValue_ = min;
                maxFlowInputValue_ = max;
            }

            /** 
             * \brief Define the continuous interval over which the flow input is sampled. 
             * @param min the minimum input values
             * @param max the maximum input values
             */
            void setJumpInputRange(std::vector<double> min, std::vector<double> max)
            {
                int size = max.size() > min.size() ? max.size() : min.size();
                if (min.size() != max.size())
                    throw Exception("Max input value (maxJumpInputValue_) and min input value (minJumpInputValue) must be of the same size");
                for (int i = 0; i < size; i++)
                {
                    if (min.at(i) > max.at(i))
                        throw Exception("Max input value must be greater than or equal to min input value");
                }
                minJumpInputValue_ = min;
                maxJumpInputValue_ = max;
            }

            /** 
             * \brief Set the maximum time allocated to a full continuous simulator step.
             * @param tM the maximum time allocated. Must be greater than 0, and greater than than the time
             * allocated to a single continuous simulator call. 
             */
            void setTm(double tM)
            {
                if (tM <= 0)
                    throw Exception("Maximum flow time per propagation step must be greater than 0");
                if (!flowStepDuration_)
                {
                    if (tM < flowStepDuration_)
                        throw Exception("Maximum flow time per propagation step must be greater than or equal to the length of time for each flow "
                                        "integration step (flowStepDuration_)");
                }
                tM_ = tM;
            }

            /** 
             * \brief Set the time allocated to a single continuous simulator call, within the full period of a continuous simulator step. 
             * @param duration the time allocated per simulator step. Must be
             * greater than 0 and less than the time allocated to a full continuous
             * simulator step. 
             */
            void setFlowStepDuration(double duration)
            {
                if (duration <= 0)
                    throw Exception("Flow step length must be greater than 0");
                if (!tM_)
                {
                    if (tM_ < duration)
                        throw Exception("Flow step length must be less than or equal to the maximum flow time per propagation step (Tm)");
                }
                flowStepDuration_ = duration;
            }

            /** 
             * \brief Set distance tolerance from goal state. 
             * @param tolerance must be greater than 0. 
             */
            void setGoalTolerance(double tolerance)
            {
                if (tolerance < 0)
                    throw Exception("Goal tolerance must be greater than or equal to 0");
                tolerance_ = tolerance;
            }

            /** 
             * \brief Define the jump set
             * @param jumpSet the jump set associated with the hybrid system. 
             */
            void setJumpSet(std::function<bool(Motion *motion)> jumpSet)
            {
                jumpSet_ = jumpSet;
            }

            /** 
             * \brief Define the flow set
             * @param flowSet the flow set associated with the hybrid system.
             */
            void setFlowSet(std::function<bool(Motion *motion)> flowSet)
            {
                flowSet_ = flowSet;
            }

            /** 
             * \brief Define the unsafe set
             * @param unsafeSet the unsafe set associated with the hybrid system.
             */
            void setUnsafeSet(std::function<bool(Motion *motion)> unsafeSet)
            {
                unsafeSet_ = unsafeSet;
            }

            /** 
             * \brief Define the distance measurement function
             * @param function the distance function associated with the motion planning problem.
             */
            void setDistanceFunction(std::function<double(base::State *, base::State *)> function)
            {
                distanceFunc_ = function;
            }

            /** 
             * \brief Define the discrete dynamics simulator
             * @param function the discrete simulator associated with the hybrid system.
             */
            void setDiscreteSimulator(std::function<base::State *(base::State *curState, const control::Control *u, base::State *newState)> function)
            {
                discreteSimulator_ = function;
            }

            /** 
             * \brief Define the continuous dynamics simulator
             * @param function the continuous simulator associated with the hybrid system.
             */
            void setContinuousSimulator(std::function<base::State *(const control::Control *u, base::State *curState, double tFlowMax,
                                                                    base::State *newState)>
                                            function)
            {
                continuousSimulator_ = function;
            }

            /** \brief Simulates the dynamics of the multicopter when in flow regime, with no nonnegligble forces other than input acceleration. */
            std::function<ompl::base::State *(const control::Control *control, ompl::base::State *x_cur, double tFlow, ompl::base::State *new_state)> continuousSimulator = [this](const control::Control *control, base::State *x_cur, double tFlow, base::State *new_state)
            {
                siC_->getStatePropagator()->propagate(x_cur, control, tFlow, new_state);
                return new_state;
            };

            /** 
             * \brief Define the collision checker
             * @param function the collision checker associated with the state space. Default is a point-by-point collision checker.
             */
            void setCollisionChecker(std::function<bool(Motion *motion, std::function<bool(Motion *motion)> obstacleSet,
                                                        base::State *newState, double *collisionTime)>
                                         function)
            {
                collisionChecker_ = function;
            }

            /** 
             * \brief Set the flow input sampling mode. 
             * @par See https://github.com/ompl/ompl/blob/main/src/ompl/util/RandomNumbers.h for details on each available mode.
             * @param mode the sampling mode.
             * @param inputs the required parameters for the given sampling mode. Not all modes require parameters.
             */
            void setFlowInputSamplingMode(inputSamplingMethods_ mode, std::vector<double> inputs)
            {
                inputSamplingMethod_ = mode;

                unsigned int targetParameterCount = 0;

                switch (mode)
                {
                case UNIFORM_INT:
                    targetParameterCount = 2;
                    getRandFlowInput_ = [this](int i)
                    { return randomSampler_->uniformInt(minFlowInputValue_[i], maxFlowInputValue_[i]); };
                    break;
                case GAUSSIAN_REAL:
                    targetParameterCount = 2;
                    getRandFlowInput_ = [&](int i)
                    { return randomSampler_->gaussian(inputs[0], inputs[1]); };
                    break;
                case HALF_NORMAL_REAL:
                    targetParameterCount = 2; // Can also be three, if want to specify focus,
                                              // which defaults to 3.0
                    getRandFlowInput_ = [&](int i)
                    { return randomSampler_->halfNormalReal(inputs[0], inputs[1], inputs[2]); };
                    break;
                default:
                    targetParameterCount = 2;
                    getRandFlowInput_ = [this](int i)
                    { return randomSampler_->uniformReal(minFlowInputValue_[i], maxFlowInputValue_[i]); };
                }

                if (inputs.size() == targetParameterCount || (mode == HALF_NORMAL_INT && targetParameterCount == 3))
                    inputSamplingParameters_ = inputs;
                else
                    throw Exception("Invalid number of input parameters for input sampling mode.");
            }

            /** 
             * \brief Set the jump input sampling mode. 
             * @par See https://github.com/ompl/ompl/blob/main/src/ompl/util/RandomNumbers.h for details on each available mode.
             * @param mode the sampling mode.
             * @param inputs the required parameters for the given sampling mode. Not all modes require parameters.
             */
            void setJumpInputSamplingMode(inputSamplingMethods_ mode, std::vector<double> inputs)
            {
                inputSamplingMethod_ = mode;

                unsigned int targetParameterCount = 0;

                switch (mode)
                {
                case UNIFORM_INT:
                    targetParameterCount = 2;
                    getRandJumpInput_ = [this](int i)
                    { return randomSampler_->uniformInt(minJumpInputValue_[i], maxJumpInputValue_[i]); };
                    break;
                case GAUSSIAN_REAL:
                    targetParameterCount = 2;
                    getRandJumpInput_ = [&](int i)
                    { return randomSampler_->gaussian(inputs[0], inputs[1]); };
                    break;
                case HALF_NORMAL_REAL:
                    targetParameterCount = 2; // Can also be three parameters, if want to specify focus,
                                              // which defaults to 3.0
                    getRandJumpInput_ = [&](int i)
                    { return randomSampler_->halfNormalReal(inputs[0], inputs[1], inputs[2]); };
                    break;
                default:
                    targetParameterCount = 2;
                    getRandJumpInput_ = [this](int i)
                    { return randomSampler_->uniformReal(minJumpInputValue_[i], maxJumpInputValue_[i]); };
                }

                if (inputs.size() == targetParameterCount || (mode == HALF_NORMAL_INT && targetParameterCount == 3))
                    inputSamplingParameters_ = inputs;
                else
                    throw Exception("Invalid number of input parameters for input sampling mode.");
            }

            /** \brief Set a different nearest neighbors datastructure */
            template <template <typename T> class NN>
            void setNearestNeighbors(void)
            {
                if (nn_ && nn_->size() != 0)
                    OMPL_WARN("Calling setNearestNeighbors will clear all states.");
                clear();
                nn_ = std::make_shared<NN<Motion *>>();
                setup();
            }

            /** \brief Check if all required parameters have been set. */
            void checkMandatoryParametersSet(void) const
            {
                if (!discreteSimulator_)
                    throw Exception("Jump map not set");
                if (!continuousSimulator_)
                    throw Exception("Flow map not set");
                if (!flowSet_)
                    throw Exception("Flow set not set");
                if (!jumpSet_)
                    throw Exception("Jump set not set");
                if (!unsafeSet_)
                    throw Exception("Unsafe set not set");
                if (!tM_)
                    throw Exception("Max flow propagation time (Tm) no set");
                // if (maxJumpInputValue_.size() == 0)
                //     throw Exception("Max input value (maxJumpInputValue) not set");
                // if (minJumpInputValue_.size() == 0)
                //     throw Exception("Min input value (minJumpInputValue) not set");
                // if (maxFlowInputValue_.size() == 0)
                //     throw Exception("Max input value (maxFlowInputValue) not set");
                // if (minFlowInputValue_.size() == 0)
                //     throw Exception("Min input value (minFlowInputValue) not set");
                // if (!flowStepDuration_)
                //     throw Exception("Flow step length (flowStepDuration_) not set");
            }

            double getFlowInput(const control::Control *control)
            {
                return control->as<CompoundControl>()->components[0]->as<RealVectorControlSpace::ControlType>()->values[0];
            }

            double getJumpInput(const control::Control *control)
            {
                return control->as<CompoundControl>()->components[1]->as<RealVectorControlSpace::ControlType>()->values[0];
            }

            void setFlowInput(control::Control *control, double value)
            {
                control->as<CompoundControl>()->components[0]->as<RealVectorControlSpace::ControlType>()->values[0] = value;
            }

            void setJumpInput(control::Control *control, double value)
            {
                control->as<CompoundControl>()->components[1]->as<RealVectorControlSpace::ControlType>()->values[0] = value;
            }

        protected:

            const static ompl::control::Control *getFlowControl(const ompl::control::Control *control)
            {
                return control->as<CompoundControl>()->as<ompl::control::Control>(0);
            }

            const static ompl::control::Control *getJumpControl(const ompl::control::Control *control)
            {
                return control->as<CompoundControl>()->as<ompl::control::Control>(1);
            }

            /** 
             * \brief Random sampler for the full vector of flow input. 
             * @return a vector of inputs (as doubles) sampled
             */
            std::function<std::vector<double>(void)> sampleFlowInputs_ = [this](void)
            {
                std::vector<double> u;
                for (unsigned int i = 0; i < maxFlowInputValue_.size(); i++)
                    u.push_back(getRandFlowInput_(i));
                return u;
            };

            /** 
             * \brief Random sampler for the full vector of jump input. 
             * @return a vector of inputs (as doubles) sampled
             */
            std::function<std::vector<double>(void)> sampleJumpInputs_ = [this](void)
            {
                std::vector<double> u;
                for (unsigned int i = 0; i < maxJumpInputValue_.size(); i++)
                    u.push_back(getRandJumpInput_(i));
                return u;
            };

            /** 
             * \brief Random sampler for one value of flow input. 
             * @return a single, double-valued sampled input
             */
            std::function<double(int i)> getRandFlowInput_ = [this](int i)
            { return randomSampler_->uniformReal(minFlowInputValue_[i], maxFlowInputValue_[i]); };

            /** 
             * \brief Random sampler for one value of jump input. 
             * @return a single, double-valued sampled input
             */
            std::function<double(int i)> getRandJumpInput_ = [this](int i)
            { return randomSampler_->uniformReal(minJumpInputValue_[i], maxJumpInputValue_[i]); };

            /// \brief The most recent goal motion.  Used for PlannerData computation
            Motion *lastGoalMotion_{nullptr};

            /** \brief Runs the initial setup tasks for the tree. */
            void initTree(void);

            /** 
             * \brief Sample the random motion. 
             * @param randomMotion The motion to be initialized
             */
            void randomSample(Motion *randomMotion);

            /// \brief A nearest-neighbors datastructure containing the tree of motions
            std::shared_ptr<NearestNeighbors<Motion *>> nn_;

            /**
             * The following are all customizeable parameters,
             * and affect how @b cHyRRT generates trajectories.
             * Customize using setter functions above. */

            /**
             * \brief Collision checker. Default is point-by-point collision checking using the jump set.
             * @param motion The motion to check for collision
             * @param obstacleSet A function that returns true if the motion's solution pair intersects with the obstacle set
             * @param ts The start time of the motion. Default is -1.0
             * @param tf The end time of the motion. Default is -1.0
             * @param newState The collision state (if a collision occurs)
             * @param collisionTime The time of collision (if a collision occurs). If no collision occurs, this value is -1.0
             * @return true if a collision occurs, false otherwise
             */
            std::function<bool(Motion *motion,
                               std::function<bool(Motion *motion)> obstacleSet, base::State *newState, double *collisionTime)>
                collisionChecker_ =
                    [this](Motion *motion,
                           std::function<bool(Motion *motion)> obstacleSet, base::State *newState, double *collisionTime = new double(-1.0)) -> bool
            {

                if (obstacleSet(motion))
                    return true;
                return false;
            };

            /// \brief Name of input sampling method, default is "uniform"
            inputSamplingMethods_ inputSamplingMethod_{UNIFORM_01};

            /// \brief Control Sampler
            control::DirectedControlSamplerPtr controlSampler_;

            /// \brief The base::SpaceInformation cast as control::SpaceInformation, for convenience
            control::SpaceInformation *siC_;

            /** 
             * \brief Compute distance between states, default is Euclidean distance 
             * @param state1 The first state
             * @param state2 The second state
             * @return The distance between the two states
             */
            std::function<double(base::State *state1, base::State *state2)> distanceFunc_ = [this](base::State *state1, base::State *state2) -> double
            {
                return si_->distance(state1, state2);
            };

            /// \brief The maximum flow time for a given flow propagation step. Must be set by the user.
            double tM_;

            /// \brief The distance tolerance from the goal state for a state to be regarded as a valid final state. Default is .1
            double tolerance_{.1};

            /// \brief The minimum step length for a given flow propagation step. Default value is 1e-6
            double minStepLength = 1e-06;

            /// \brief The flow time for a given integration step, within a flow propagation step. Must be set by user.
            double flowStepDuration_;

            /// \brief Minimum flow input values
            std::vector<double> minFlowInputValue_;

            /// \brief Maximum flow input values
            std::vector<double> maxFlowInputValue_;

            /// \brief Minimum jump input values
            std::vector<double> minJumpInputValue_;

            /// \brief Maximum jump input values
            std::vector<double> maxJumpInputValue_;

            /** 
             * \brief Simulator for propagation under jump regime
             * @param curState The current state
             * @param u The input
             * @param newState The newly propagated state
             * @return The newly propagated state
             */
            std::function<base::State *(base::State *curState, const control::Control *u, base::State *newState)> discreteSimulator_;

            /** 
             * \brief Function that returns true if a motion intersects with the jump set, and false if not. 
             * @param motion The motion to check
             * @return True if the state is in the jump set, false if not
             */
            std::function<bool(Motion *motion)> jumpSet_;

            /** 
             * \brief Function that returns true if a motion intersects with the flow set, and false if not. 
             * @param motion The motion to check
             * @return True if the state is in the flow set, false if not
             */
            std::function<bool(Motion *motion)> flowSet_;

            /** 
             * \brief Function that returns true if a motion intersects with the unsafe set, and false if not. 
             * @param motion The motion to check
             * @return True if the state is in the unsafe set, false if not
             */
            std::function<bool(Motion *motion)> unsafeSet_;

            /** 
             * \brief Simulator for propagation under flow regime
             * @param input The input
             * @param curState The current state
             * @param tFlowMax The random maximum flow time
             * @param newState The newly propagated state
             * @return The newly propagated state
             */
            std::function<base::State *(const control::Control *u, base::State *curState, double tFlowMax, base::State *newState)> continuousSimulator_;

            /// \brief Random sampler for the input. Default constructor always seeds a different value, and returns a uniform real distribution.
            RNG *randomSampler_ = new RNG();

            /** 
             * \brief Construct the path, starting at the last edge. 
             * @param lastMotion The last motion in the solution
             * @return the planner status (APPROXIMATE, EXACT, or UNKNOWN)
             */
            base::PlannerStatus constructSolution(Motion *lastMotion);

            /// \brief Input sampling parameters
            std::vector<double> inputSamplingParameters_{};

            /// \brief State sampler
            base::StateSamplerPtr sampler_;
        };
    }
}

#endif