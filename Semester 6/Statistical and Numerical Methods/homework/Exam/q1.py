import numpy as np
import random
import matplotlib.pyplot as plt

class Building:
    def __init__(self, floors, elevator_rest_place=None):
        self.floors = floors
        self.cost_per_floor = 1
        self.elevator_rest_place = elevator_rest_place
    
    def get_travel_cost(self, start_floor, end_floor):
        if start_floor < 0 or start_floor > self.floors or end_floor < 0 or end_floor > self.floors:
            raise ValueError("Invalid floor number")
        return (abs(self.elevator_rest_place - start_floor) + abs(end_floor - start_floor) + abs(self.elevator_rest_place - end_floor)) * self.cost_per_floor

def simulate_day_one_person(building, travel_floors=None):
    if travel_floors is None:
        raise ValueError("travel_floors must be provided for the person")
    total_cost = 0
    current_floor = 0
    # travel to all the floors
    for destination in travel_floors:
        total_cost += building.get_travel_cost(current_floor, destination)
        current_floor = destination
    return total_cost

def simulate_day(building, simulation_people=5000):
    # simulate 5000 people
    total_cost = 0
    for _ in range(simulation_people):
        travel_floors = []
        current = 0
        while len(travel_floors) < 4:
            next_floor = random.randint(1, building.floors)
            # only ensure the next floor is not the current (last) floor
            if next_floor != current:
                travel_floors.append(next_floor)
                current = next_floor
        travel_floors.append(0)  # ensure the last floor is ground floor
        total_cost += simulate_day_one_person(building, travel_floors)

    return total_cost

def main():
    floors = 10
    simulation_runs = 100
    simulation_people = 5000
    elevator_rest_places_candidates = np.linspace(0, floors, num=11)
    rest_costs = {place: 0 for place in elevator_rest_places_candidates}
    for _ in range(simulation_runs):
        for elevator_rest_place in elevator_rest_places_candidates:
            building = Building(floors, elevator_rest_place)
            total_cost = simulate_day(building, simulation_people)
            rest_costs[elevator_rest_place] += total_cost
    average_costs = {place: cost / simulation_runs / simulation_people for place, cost in rest_costs.items()}
    optimal_rest_place = min(average_costs.items(), key=lambda kv: kv[1])[0]
    print(f"Optimal elevator rest place: {optimal_rest_place}, Minimum cost: {average_costs[optimal_rest_place]}")
    plt.plot(list(average_costs.keys()), list(average_costs.values()))
    plt.xlabel("Elevator Rest Place")
    plt.ylabel("Total Cost")
    plt.title("Cost vs Elevator Rest Place")
    plt.show()

if __name__ == "__main__":
    main()