#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <queue>
#include <chrono>
#include <fstream>


// All lengths in meters, times in seconds

class Trolley {

public:

    Trolley(){}

    Trolley(int capacity, double max_velocity, double acceleration,
            double length, double start_coord, double path_length) {
        this->capacity = capacity;
        this->num_people = 0;

        this->max_velocity = max_velocity;
        this->velocity = 0;
        this->acceleration = acceleration;
        this->coordinate = start_coord;

        this->length = length;
        this->path_length = path_length;

        stop_index = -1;
    }

    void step(double dt, double max_shift){
        if (velocity + acceleration * dt < max_velocity && max_shift != 0)
            velocity += acceleration * dt;

        coordinate += std::min(velocity * dt, max_shift);

        if (coordinate > path_length) coordinate -= path_length;
    }

    void exitingPeople (std::uniform_real_distribution<double>& distribution, std::default_random_engine& engine) {
        int people_to_exit = int(num_people * distribution(engine));

        num_people -= people_to_exit;
    }

    int loadingPeople(int num_people_come) {
        int rest = num_people + num_people_come - capacity;

        if (rest < 0) {
            num_people += num_people_come;
            return 0;
        }

        num_people = capacity;
        return rest;

    }

    double getCoord() const {
        return coordinate;
    }

    double getLength() const {
        return length;
    }

    double getNumPeople() const {
        return num_people;
    }

    int getCapacity() const {
        return capacity;
    }

    int getStopIndex() const {
        return stop_index;
    }

    void onStop(int index) {
        stop_index = index;
    }

    void leaveStop() {
        stop_index = -1;
    }

private:
    int capacity;
    int num_people;

    double max_velocity;
    double velocity;
    double acceleration;
    double coordinate;

    double length;
    double path_length;

    int stop_index;
};

bool operator< (const Trolley& lhs, const Trolley& rhs) {
    return lhs.getCoord() < rhs.getCoord();
}



class Stop {

public:

    Stop()= default;

    Stop(double coordinate) {
        this->coordinate = coordinate;

        num_people = 0;
    }

    void peopleCome(std::bernoulli_distribution& distribution, std::default_random_engine& engine) {
        num_people += int(distribution(engine));
    }

    void peopleLeft(int number) {
        num_people -= number;
    }

    int getNumPeople() const {
        return num_people;
    }

    double getCoord() const {
        return coordinate;
    }

private:

    int num_people{};

    double coordinate{};
};


class Controller {

public:

    Controller()= default;

    Controller(double path_length, const std::vector<double>& trolleys_coords, const std::vector<double>& stops_coords,
               int trolley_capacity, double trolley_max_velocity,
               double trolley_acceleration, double trolley_length, double step_time, double chance_person_per_step) {

        this->path_length = path_length;

        for (auto x : trolleys_coords) {
            trolleys.emplace_back(trolley_capacity, trolley_max_velocity,
                                  trolley_acceleration, trolley_length, x,path_length);
        }

        std::sort(std::begin(trolleys), std::end(trolleys));

        for (auto x : stops_coords) {
            stops.emplace_back(x);
        }

        this->step_time = step_time;

        min_distances.resize(trolleys_coords.size());

        engine = std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());

        distribution = std::bernoulli_distribution(chance_person_per_step);

        exiting_dist = std::uniform_real_distribution<double>(0, 0.8);

        out.open("data.txt");

    }

    void run(int number_steps) {

        int step_counter = 0;

        while (step_counter < number_steps) {
            step_counter++;
            movement();
            logic();
            if (step_counter % 60 == 0) {
                print_results(step_counter);
                collect_data(int(step_counter / 60));
            }
        }
    }

    void logic() {

        for(auto& stop : stops) {
            stop.peopleCome(distribution, engine);
        }

        for (int i = 0; i < trolleys.size(); i++) {
            double min_distance = 0;
            if (i < trolleys.size() - 1) min_distance = dist(trolleys[i].getCoord(),
                                       trolleys[i + 1].getCoord() - trolleys[i + 1].getLength());

            else min_distance = dist(trolleys[i].getCoord(),
                                     trolleys[0].getCoord() - trolleys[0].getLength());

            if (trolleys[i].getStopIndex() != -1) {
                if (!processStops(i)) {
                    trolleys[i].leaveStop();
                }
                else {
                    min_distance = 0;
                }
            }
            else {
                for (int j = 0; j < stops.size(); j++) {

                    if (stops[j].getCoord() > trolleys[i].getCoord() - 1 &&
                            stops[j].getCoord() < trolleys[i].getCoord() + 1) {
                        trolleys[i].onStop(j);
                        trolleys[i].exitingPeople(exiting_dist, engine);
                        min_distance = 0;
                    }
                    else {
                        min_distance = std::min(min_distance,
                                                dist(trolleys[i].getCoord(), stops[j].getCoord()));
                    }

                }
            }

            min_distances[i] = min_distance;
        }

    }

    void movement() {

        for (int i = 0; i < trolleys.size(); i++) {
            trolleys[i].step(step_time, min_distances[i]);
        }
        std::sort(std::begin(trolleys), std::end(trolleys));
    }

    double dist(double coord_1, double coord_2) const {
        if (coord_1 > path_length) coord_1 -= path_length;
        if (coord_2 > path_length) coord_2 -= path_length;

        double distance = std::abs(coord_1 - coord_2);
        return std::min(distance, path_length - distance);
    }

    bool processStops(int trolley_index) {

        int stop_index = trolleys[trolley_index].getStopIndex();

        if (stops[stop_index].getNumPeople() == 0 ||
            trolleys[trolley_index].getNumPeople() == trolleys[trolley_index].getCapacity()) return false;

        stops[stop_index].peopleLeft(people_per_step_time);
        trolleys[trolley_index].loadingPeople(people_per_step_time);

        return true;
    }

    void print_results(int step) {
        std::cout << "time: " << step * step_time << "\n";

        for (auto& trolley : trolleys) {
            std::cout << "coordinate: " << trolley.getCoord() << ", number people: " << trolley.getNumPeople() << "\n";
        }

        for(auto& stop : stops) {
            std::cout << "number people on stop: " << stop.getNumPeople() << "\n";
        }

        std::cout << "\n\n";
    }

    void collect_data(int minutes) {
        out << minutes << " " << calculate_average_min_dist() << std::endl;
    }

    double calculate_average_min_dist() {
        double sum = std::min(dist(trolleys[0].getCoord(), trolleys[trolleys.size() - 1].getCoord()),
                              dist(trolleys[0].getCoord(), trolleys[1].getCoord()));

        sum += std::min(dist(trolleys[0].getCoord(), trolleys[trolleys.size() - 1].getCoord()),
                        dist(trolleys[trolleys.size() - 2].getCoord(),
                             trolleys[trolleys.size() - 1].getCoord()));

        for (int i = 1; i < trolleys.size() - 1; i++) {
            sum += std::min(dist(trolleys[i].getCoord(), trolleys[i - 1].getCoord()),
                            dist(trolleys[i].getCoord(), trolleys[i + 1].getCoord()));
        }

        return sum / trolleys.size();
    }


private:
    double path_length{};
    double step_time{};

    const int people_per_step_time = 1;

    std::vector<Trolley> trolleys;
    std::vector<Stop> stops;
    std::vector<double> min_distances;


    std::bernoulli_distribution distribution;
    std::default_random_engine engine;
    std::uniform_real_distribution<double> exiting_dist;

    std::ofstream out;
};

int main() {
    double path_length = 8000;
    double trolleys_length = 20;
    int trolleys_capacity = 30;
    double max_vel = 20;
    double acceleration = 0.4;
    double chance_spawn_person_per_step_time = 0.13;
    double step_time = 1;

    std::vector<double> trolley_coords = {50, 1000, 2000, 3000, 4000, 5000, 6000, 7000};

    std::vector<double> stops_coords = {500, 1500, 2500, 3500, 4500, 5500, 6500, 7500};


    Controller controller(path_length, trolley_coords, stops_coords, trolleys_capacity, max_vel,
                          acceleration, trolleys_length, step_time, chance_spawn_person_per_step_time);


    int number_steps = 100000;

    controller.run(number_steps);

    return 0;
}
