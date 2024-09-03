// -----------------------------------------------------------------------------
//
// Copyright (C) 2024 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
//
// A simulation of a 2D wound healing assay with walls of cells
//

#ifndef SCRATCH_WOUND_ASSAY_H_   // To prevents multiple inclusions of the same header file
#define SCRATCH_WOUND_ASSAY_H_

#include "biodynamo.h"  // To include BioDynaMo framework which provides core functionalities

namespace bdm {  // To create a namespace "bdm" to encapsulate the code and avoid conflicts

// MyCell class defines a custom cell type extending the existing "Cell" class
class MyCell : public Cell {
  BDM_AGENT_HEADER(MyCell, Cell, 1);  // Macro to handle boilerplate for agents (cells in this case)

 public:
  MyCell() {}  // This is the default constructor
  explicit MyCell(const Real3& position) : Base(position) {}  // Constructor that takes a position
  virtual ~MyCell() {}  // Destructor

  // To initialize function called when a new cell is created via division
  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);  // To call parent class's initialization
    if (auto* mother = dynamic_cast<MyCell*>(event.existing_agent)) {
      // If the new cell is a daughter from a division, it inherits properties from the mother cell
      can_divide_ = mother->can_divide_;
      initial_side_ = mother->initial_side_;
      migrating_ = mother->migrating_;
    }
  }

  // Setter and getter for the ability to divide
  void SetCanDivide(bool d) { can_divide_ = d; }
  bool GetCanDivide() const { return can_divide_; }

  // Setter and getter for initial side (left or right surface in the wound model)
  void SetInitialSide(int side) { initial_side_ = side; }
  int GetInitialSide() const { return initial_side_; }

  // Setter and getter for the migrating status
  void SetMigrating(bool migrating) { migrating_ = migrating; }
  bool IsMigrating() const { return migrating_; }

 private:
  bool can_divide_ = true;  // Determines if the cell can divide
  int initial_side_;  // Identifies if the cell is on the left (0) or right (1) side of the wound
  bool migrating_ = false;  // Flag to indicate if the cell is migrating into the wound gap
};

// WoundHealingBehavior class defines how cells behave during wound healing
struct WoundHealingBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(WoundHealingBehavior, Behavior, 1);  // Macro for agent behaviors

  WoundHealingBehavior() { AlwaysCopyToNew(); }  // Ensure behavior is copied to new agents (cells)
  virtual ~WoundHealingBehavior() {}

  // Run function defines the behavior that is executed on each simulation step
  void Run(Agent* agent) override {
    if (auto* cell = dynamic_cast<MyCell*>(agent)) {  // To ensure the agent is a MyCell type
      auto* sim = Simulation::GetActive();  // To access the current simulation instance
      auto* ctxt = sim->GetExecutionContext();  // To get execution context to manage agent interactions

      bool gap_detected = false;  // This will track if there's a gap in the cell layer
      int side = cell->GetInitialSide();  // To get the side the cell is on (left or right)

      // We assume cell is on the edge initially (until we find neighbors)
      bool is_edge_cell = true;

      // Lambda function to check for neighboring cells within a certain set range along the x-axis
      // The range can be changed according to need
      auto check_neighbors = L2F([&](Agent* neighbor, real_t squared_distance) {
        if (side == 0) {  // Left side cells look for neighbors to their right
          if (neighbor->GetPosition()[0] > cell->GetPosition()[0] &&
              neighbor->GetPosition()[0] < cell->GetPosition()[0] + 20) {
            is_edge_cell = false;  // If neighbor found on the right, it's not an edge cell
          }
        } else if (side == 1) {  // Right side cells look for neighbors to their left
          if (neighbor->GetPosition()[0] < cell->GetPosition()[0] &&
              neighbor->GetPosition()[0] > cell->GetPosition()[0] - 20) {
            is_edge_cell = false;  // If neighbor found on the left, it's not an edge cell
          }
        }
      });

      // -----------------------------------------------------------------------------
      // Cells on the left side will look for neighbors on the right and cells on the right side will on look left
      // This prevents cells from considering the spaces of the simualation boundary that was not used which created so the size of the agent surface
      // can be adjusted, is not considered a wound 
      // This prevents the cells from thinking the empty space to above, below and to the left of the left surface and right of the right surface are wounds
      // -----------------------------------------------------------------------------



      // To execute the search for neighboring cells within a radius of 30 units
      ctxt->ForEachNeighbor(check_neighbors, *cell, 30);

      // If the cell is on the edge, we check for a gap in the wound area
      if (is_edge_cell) {
        auto search_gap = L2F([&](Agent* neighbor, real_t squared_distance) {
          if (side == 0) {  // Left side checks for gaps to the right
            if (neighbor->GetPosition()[0] > cell->GetPosition()[0] &&
                neighbor->GetPosition()[0] < cell->GetPosition()[0] + 30) {
              gap_detected = true;  // Neighbor detected, no gap
            }
          } else if (side == 1) {  // Right side checks for gaps to the left
            if (neighbor->GetPosition()[0] < cell->GetPosition()[0] &&
                neighbor->GetPosition()[0] > cell->GetPosition()[0] - 30) {
              gap_detected = true;  // Neighbor detected, no gap
            }
          }
        });

        // Search again within a radius of 40 units for a gap
        ctxt->ForEachNeighbor(search_gap, *cell, 40);

        // If no gap is found, and the cell can divide, it creates a daughter cell
        if (!gap_detected && cell->GetCanDivide()) {
          auto daughter_position = cell->GetPosition();  // Get the current position of the cell

          // To position the daughter cell closer to the gap
          if (side == 0) {  // If left side, place the daughter to the right
            daughter_position[0] += cell->GetDiameter();
          } else if (side == 1) {  // If right side, place the daughter to the left
            daughter_position[0] -= cell->GetDiameter();
          }

          // To create the new daughter cell and position it
          auto* daughter = dynamic_cast<MyCell*>(cell->Divide());
          daughter->SetPosition(daughter_position);  // Set the position of the daughter cell

          // To mark the daughter as migrating towards the gap
          daughter->SetMigrating(true);
        }
      }

      // Behavior for cells marked as migrating
      if (cell->IsMigrating()) {
        // Cells grow until they reach a certain diameter (size)
        // 8.0 is chosen arbitrarilly 
        if (cell->GetDiameter() < 8.0) {
          cell->ChangeVolume(400);  // To increase the cell's volume (size)
        } else {
          // Once fully grown, stop migration
          cell->SetMigrating(false);
        }
      }
    }
  }
};

// Main simulation function
inline int Simulate(int argc, const char** argv) {
  // To setup the simulation parameters (space boundaries)
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;  // To enclose the simulation space
    param->min_bound = -1000;  // To set the minimum space boundary
    param->max_bound = 1000;  // To set the maximum space boundary
  }; 

  // To create the simulation with parameters
  Simulation simulation(argc, argv, set_param);
  auto* ctxt = simulation.GetExecutionContext();  // To get the execution context

  // Function to create a new cell at a specific position and side (left/right)
  auto create_agent = [](const Real3& position, int side) {
    MyCell* cell = new MyCell(position);  // To create a new cell
    cell->SetDiameter(7.5);  // To set the initial cell diameter
    cell->AddBehavior(new WoundHealingBehavior());  // To assign the wound healing behavior
    cell->SetInitialSide(side);  // To set the cell's side (left or right)
    return cell;
  };

  // To define a flat surface at z = 0
  auto f = [](const real_t* x, const real_t* params) {
    return 0.0;
  };

  // To create cells on two surfaces 
  // The space (gap) between them acts as the wound that the daughter cells need to proliferate into to heal
  // The x parameters can be changed to increase the size of the wound  
  ModelInitializer::CreateAgentsOnSurface(f, {}, -200, -150, 10, 0, 500, 10,
                                          [&](const Real3& pos) { return create_agent(pos, 0); });
  ModelInitializer::CreateAgentsOnSurface(f, {}, 100, 150, 10, 0, 500, 10,
                                          [&](const Real3& pos) { return create_agent(pos, 1); });

  // To print message that the simulation has started
  // This is purely for error handling purposes
  std::cout << "Wound healing simulation started successfully!" << std::endl;
  
  // To run the simulation for 500 timesteps
  // The timestep can be varied to nneed
  simulation.GetScheduler()->Simulate(500);

  // Print message that the simulation has completed
  std::cout << "Wound healing simulation completed successfully!" << std::endl;
  return 0;  // To exit the simulation
}

}  // namespace bdm

#endif  // SCRATCH_WOUND_ASSAY_H_
