#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

struct Interval {
    int lower, upper;

    Interval(int l, int u) : lower(l), upper(u) {}

    bool operator<(const Interval& rhs) const {
        return (lower < rhs.lower) || (lower == rhs.lower && upper < rhs.upper);
    }

    bool overlaps(const Interval& rhs) const {
        return !(upper < rhs.lower || lower > rhs.upper);
    }

    Interval merge(const Interval& rhs) const {
        return Interval(std::min(lower, rhs.lower), std::max(upper, rhs.upper));
    }
};

struct RV {
    std::string name;
    std::set<Interval> possibleValues;

    RV(const std::string& name, const std::set<Interval>& values) : name(name), possibleValues(values) {}

    bool isCompatible(const std::set<Interval>& intervals) const {
        // Check if intervals are within the RV's possible values
        for (const auto& interval : intervals) {
            bool found = false;
            for (const auto& possible : possibleValues) {
                if (interval.overlaps(possible)) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }
};

struct Event_expr;

struct Event {
    std::string name;
    RV* RVName;
    double Probability;
    std::set<Interval> Values;

    Event(const std::string& name, RV* rvName, const std::set<Interval>& intervals)
        : name(name), RVName(rvName), Values(intervals) {
        // Probability = calculateProbability(); // Uncomment and implement as needed
    }

    Event_expr operator^(const Event& rhs) const;
    Event_expr operator+(const Event& rhs) const;
    Event_expr operator|(const Event& rhs) const;
    void operator=(const Event_expr& rhs);
};

struct Event_expr {
    enum class Op { NONE, UNION, INTERSECT, CONDITIONAL };

    Op op;
    std::vector<Event> events;

    Event_expr(const Event& e) : op(Op::NONE) {
        events.push_back(e);
    }
};

// Helper function to merge intervals
std::set<Interval> mergeIntervals(const std::set<Interval>& a, const std::set<Interval>& b) {
    std::set<Interval> result = a;
    for (const auto& interval : b) {
        bool merged = false;
        for (auto it = result.begin(); it != result.end(); ) {
            if (it->overlaps(interval)) {
                Interval mergedInterval = it->merge(interval);
                it = result.erase(it);
                result.insert(mergedInterval);
                merged = true;
                break;
            } else {
                ++it;
            }
        }
        if (!merged) {
            result.insert(interval);
        }
    }
    return result;
}

// Union operator: merge intervals and check RV compatibility
Event_expr Event::operator+(const Event& rhs) const {
    if (RVName->name != rhs.RVName->name) {
        throw std::invalid_argument("Cannot union events with different random variables.");
    }

    std::set<Interval> mergedIntervals = mergeIntervals(Values, rhs.Values);

    if (!RVName->isCompatible(mergedIntervals)) {
        throw std::invalid_argument("Union of intervals is not compatible with the random variable.");
    }

    Event result("Union(" + name + ", " + rhs.name + ")", RVName, mergedIntervals);
    return Event_expr(result);
}

// Intersection operator: find overlapping intervals
Event_expr Event::operator^(const Event& rhs) const {
    if (RVName->name != rhs.RVName->name) {
        throw std::invalid_argument("Cannot intersect events with different random variables.");
    }

    std::set<Interval> intersectedIntervals;
    for (const auto& intervalA : Values) {
        for (const auto& intervalB : rhs.Values) {
            if (intervalA.overlaps(intervalB)) {
                intersectedIntervals.insert(intervalA.merge(intervalB));
            }
        }
    }

    if (!RVName->isCompatible(intersectedIntervals)) {
        throw std::invalid_argument("Intersection of intervals is not compatible with the random variable.");
    }

    Event result("Intersection(" + name + ", " + rhs.name + ")", RVName, intersectedIntervals);
    return Event_expr(result);
}

// Conditional operator: this could be defined differently depending on semantics
Event_expr Event::operator|(const Event& rhs) const {
    if (RVName->name != rhs.RVName->name) {
        throw std::invalid_argument("Cannot condition on events with different random variables.");
    }

    // Here you can define custom behavior for the conditional, for example, combining intervals.
    std::set<Interval> conditionalIntervals = mergeIntervals(Values, rhs.Values);

    if (!RVName->isCompatible(conditionalIntervals)) {
        throw std::invalid_argument("Conditional intervals are not compatible with the random variable.");
    }

    Event result("Conditional(" + name + ", " + rhs.name + ")", RVName, conditionalIntervals);
    return Event_expr(result);
}

// Assignment operator for Event_expr
void Event::operator=(const Event_expr& rhs) {
    if (rhs.events.empty()) {
        throw std::invalid_argument("Event_expr has no events to assign.");
    }

    const Event& assignedEvent = rhs.events[0];
    if (RVName->name != assignedEvent.RVName->name) {
        throw std::invalid_argument("Cannot assign Event_expr with a different random variable.");
    }

    name = assignedEvent.name;
    Values = assignedEvent.Values;
}
