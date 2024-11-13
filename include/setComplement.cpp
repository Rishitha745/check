#include <iostream>
#include <set>
#include <utility>
using namespace std;


set<pair<int,int>> complement(set<pair<int,int>>& a, set<pair<int,int>>& b) {
    set<pair<int,int>> result;
    if (a.empty()) return result;
    if (b.empty()) {
        result = a;
        return result;
    }
    auto it1 = a.begin();
    auto it2 = b.begin();
    while (it1 != a.end()) {
        if (it2 == b.end() || it1->second < it2->first) {
            result.insert(*it1);
            ++it1;
            continue;
        }
        if (it1->first > it2->second) {
            ++it2;
            continue;
        }
        if (it1->first < it2->first) {
            result.insert({it1->first, it2->first});
        }
        if (it1->second > it2->second) {
            int start = it2->second;
            int end = it1->second;
            ++it2;
            while (it2 != b.end() && it2->first <= end) {
                if (start < it2->first) {
                    result.insert({start, it2->first});
                }
                start = max(start, it2->second);
                if (start >= end) break;
                ++it2;
            }
            if (start < end) {
                result.insert({start, end});
            }
        } else {
            ++it1;
            continue;
        }
        ++it1;
    }
    
    return result;
}



set<pair<int,int>> complement2(set<pair<int,int>>& a, set<pair<int,int>>& b) {
    set<pair<int,int>> result;
    if (a.empty()) return result;
    if (b.empty()) {
        result = a;
        return result;
    }
    auto it1 = a.begin();
    auto it2 = b.begin();
    while (it1 != a.end()) {
        if (it2 == b.end() || it1->second < it2->first) {
            result.insert(*it1);
            ++it1;
            continue;
        }
        if (it1->first > it2->second) {
            ++it2;
            continue;
        }
        if (it1->first < it2->first) {
            result.insert({it1->first, it2->first});
        }
        if (it1->second > it2->second) {
            int start = it2->second;
            int end = it1->second;
            ++it2;
            while (it2 != b.end() && it2->first <= end) {
                if (start < it2->first) {
                    result.insert({start, it2->first});
                }
                start = max(start, it2->second);
                if (start >= end) break;
                ++it2;
            }
            if (start < end) {
                result.insert({start, end});
            }
        } else {
            if (it1->second == it2->second) {
                ++it2;
            }
            ++it1;
            continue;
        }
        ++it1;
    }
    
    return result;
}


int main(){
    set<pair<int,int>> a,b;
    int n=0,m=0;
    int p=0,q=0;
    cin>>n>>m;
    pair<int,int> z ;
    for ( int i=0;i<n;i++){
        cin>>p;
        cin>>q;
        z = {p,q};
        a.insert(z); 
    }
    for ( int i=0;i<m;i++){
        cin>>p;
        cin>>q;
        z = {p,q};
        b.insert(z); 
    }
    set<pair<int,int>> c = complement(a,b);
    for(auto [a,b] : c){
        cout<<a<<" "<<b<<endl;
    }

}