#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <vector>
using namespace std;
#define all(x) x.begin(),x.end()
#define ppoint pair<double, double>
#define ppoint3D pair<ppoint, double>
#define INF 1000000000

ppoint rotate_around_origin(ppoint pt, double angle){
  double s = sin(angle);
  double c = cos(angle);
  return make_pair(pt.first * c - pt.second * s, pt.first * s + pt.second * c);
}

vector<double> find_smallest_rectangle(vector<ppoint> pts){
  int N_ROTATIONS = 1000;
  double PI = 3.14159265;
  double min_area = 1e9;
  double best_angle = -1;
  ppoint pt1, pt2;
  for(int i = 0; i < N_ROTATIONS; i++){
    double angle = (2.0 * PI / N_ROTATIONS) * i;
    double maxX = -INF, maxY = -INF;
    double minX = INF, minY = INF;
    for(auto it = pts.begin(); it != pts.end(); ++it){
      ppoint tmp = rotate_around_origin(*it, angle);
      maxX = max(maxX, tmp.first);
      minX = min(minX, tmp.first);
      maxY = max(maxY, tmp.second);
      minY = min(minY, tmp.second);
    }
    double area = (maxX - minX) * (maxY - minY);
    if(area < min_area){
      min_area = area;
      best_angle = angle;
      pt1 = rotate_around_origin(make_pair(minX, maxY), -angle);
      pt2 = rotate_around_origin(make_pair(maxX, minY), -angle);
    }
  }
  assert(best_angle != -1);
  vector<double> res {pt1.first, pt1.second, pt2.first, pt2.second, best_angle};
  return res;
}

double area(ppoint a, ppoint b, ppoint c){
  return (b.first - a.first) * (c.second - a.second) - (b.second - a.second) * (c.first - a.first);
}

vector<ppoint> find_convex_hull(vector<ppoint> P){
  sort(all(P));
  vector<ppoint> upper, lower;
  upper.push_back(P[0]);
  lower.push_back(P[0]);
  for(int i = 0; i < P.size(); i++){
    if(i == P.size() - 1 || area(P[0], P[i], P.back()) > 0){
      while(upper.size() >= 2 && area(upper[upper.size()-2], upper.back(), P[i]) < 0)
        upper.pop_back();
      upper.push_back(P[i]);
    }
    if(i == P.size() - 1 || area(P[0], P[i], P.back()) < 0){
      while(lower.size() >= 2 && area(lower[lower.size()-2], lower.back(), P[i]) > 0)
        lower.pop_back();
      lower.push_back(P[i]);
    }
  }
  for(int i = lower.size() - 2; i >= 1; i--)
    upper.push_back(lower[i]);
  return upper;
}

vector<double> remove_surface_points(vector<ppoint3D> V){
  // her surface point için "point[label] = eliminated;"
  map<int, int> mp;
  for(auto it = V.begin(); it != V.end(); ++it)
      mp[int(it->first.second * 100)]++;
  int surfaceZ = 0, maxi = 0;
  for(auto it = mp.begin(); it != mp.end(); ++it)
    if(it->second > maxi){
      maxi = it->second;
      surfaceZ = it->first;
    }
    double r = 2.0;
    int s1 = 0, s2 = 0, s3 = 0;
    double a=0, b=0, c=0;
    for(auto it = V.begin(); it != V.end(); ++it){
        int x = int(it->first.second * 100);
        if(surfaceZ - r <= x && x <= surfaceZ + r) {
          a += it->first.first; b += it->first.second; c += it->second;
          s3++;
        }
        else {
          s1++;
          if(x < surfaceZ - r)
            s2++;
        }
      }

    if(s2 * 3 < s1)
      for(auto it = V.begin(); it != V.end(); ++it)
          if(int(it->first.second * 100) <= surfaceZ + r){
            a += it->first.first; b += it->first.second; c += it->second;
          }
    vector<double> res {a, b, c, double(s1), double(s2), double(s3)};
    return res;
}

int main(){

  freopen("input.txt", "r", stdin);
  int N;
  double a,b,c;
  vector<ppoint> pts;
  vector<ppoint3D> pts3D;
  cin >> N;
  while(N--){
    cin >> a >> b >> c;
    pts.push_back(make_pair(a,b));
    pts3D.push_back(make_pair(make_pair(a,b),c));
  }

  vector<double> res = find_smallest_rectangle(pts);
  printf("%.5f %.5f %.5f %.5f %.5f\n\n", res[0], res[1], res[2], res[3], res[4]);

  vector<ppoint> res2 = find_convex_hull(pts);
  for(auto it = res2.begin(); it != res2.end(); ++it)
    printf("%.5f %.5f\n", it->first, it->second);
  cout << endl;

  vector<double> res3 = remove_surface_points(pts3D);
  for(auto it = res3.begin(); it != res3.end(); ++it)
    printf("%.5f ", *it);
  cout << endl;
  return 0;
}
