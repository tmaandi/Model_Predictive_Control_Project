#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

/* Calibration flag to tune cost weightages */
#define CALIB 0

using json = nlohmann::json;

/* For converting back and forth between radians and degrees. */
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

/*
 * Checks if the SocketIO event has JSON data.
 * If there is data the JSON object in string format will be returned,
 * else the empty string "" will be returned.
 */

string hasData(string s)
{
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos)
  {
    return "";
  }
  else if (b1 != string::npos && b2 != string::npos)
  {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

/* Evaluate a polynomial. */
double polyeval(Eigen::VectorXd coeffs, double x)
{
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++)
  {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

/*
 * Fit a polynomial.
 * Adapted from
 * https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
 */

Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order)
{
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++)
  {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++)
  {
    for (int i = 0; i < order; i++)
    {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main(int argc, const char *argv[])
{
  uWS::Hub h;

  if (CALIB == 1)
  {
   cte_w = stod(argv[1]);
   psie_w = stod(argv[2]);
   delta_w = stod(argv[3]);
  }
  else
  {
    cte_w = 2;
    psie_w = 100;
    delta_w = 5;
  }

  /* MPC is initialized here! */
  MPC mpc;

  h.onMessage([&mpc, &Lf](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode)
  {

    latency = 0.1; /* 100 ms */

    /*
     * "42" at the start of the message means there's a websocket message event.
     * The 4 signifies a websocket message
     * The 2 signifies a websocket event
     */
    string sdata = string(data).substr(0, length);
    std::cout << sdata << std::endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2')
    {
      string s = hasData(sdata);
      if (s != "")
      {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          std::vector<double> ptsx = j[1]["ptsx"];
          std::vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double delta = j[1]["steering_angle"];
          double throttle = j[1]["throttle"];

          double acceleration = 0.5 * throttle; /* Assuming lag in throttle to resulting accl. */

          /*
           * Finding relative diff. and trasforming from
           * global to vehicle coordinates
           */
          for (unsigned int i = 0; i < ptsx.size(); ++i)
          {
            double loc_x = ptsx[i] - px;
            double loc_y = ptsy[i] - py;

            /* Coordinate Transformation */
            ptsx[i] = loc_x * cos(psi) + loc_y * sin (psi);
            ptsy[i] =  - loc_x * sin(psi) + loc_y * cos (psi);
          }

          /* Update to vehicle frame of reference */
          px = 0;
          py = 0;
          psi = 0;

          /* Data type conversion */
          Eigen::VectorXd ptx = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptsx.data(), ptsx.size());
          Eigen::VectorXd pty = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptsy.data(), ptsy.size());

          /* Fit a polynomial to the waypoint x and y coordinates */
          auto coeffs = polyfit(ptx, pty, 3);
          /* Calculate the cross track error */
          double cte = (py - polyeval(coeffs, px));
          /* Calculate the orientation error */
          double epsi = (psi - atan(coeffs[1] + 2 * coeffs[2] * px + 3 * coeffs[3] * px * px));

          /* State Vector */
          Eigen::VectorXd state(6);

          /*
           * Correcting for sensor latency
           */
          px = px + v * cos(psi) * latency;
          py = py + v * sin(psi) * latency;
          psi = psi - v * delta/Lf * latency;
          v = v + acceleration * latency;

          state << px, py, psi, v, cte, epsi;

          vector<double> optim_output;
          double steer_value;
          double throttle_value;

          optim_output = mpc.Solve(state, coeffs);

          steer_value = optim_output[0];
          throttle_value = optim_output[1];

          json msgJson;
          msgJson["steering_angle"] = -steer_value;
          msgJson["throttle"] = throttle_value;

          /* Display the MPC predicted trajectory */
          vector<double> mpc_x_vals = {optim_output[2], optim_output[3], optim_output[4],
              optim_output[5], optim_output[6],optim_output[7]};
          vector<double> mpc_y_vals = {optim_output[8], optim_output[9], optim_output[10],
              optim_output[11], optim_output[12],optim_output[13]};

          /*
           * Here, points are in reference to the vehicle's coordinate system
           * the points in the simulator are connected by a Green line
           */

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;


          /* Display the waypoints/reference line */
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (int i = 0; i < 50; ++i)
          {
            next_y_vals.push_back(polyeval(coeffs, i));
            next_x_vals.push_back(i);
          }

          /*
           * Way points are in reference to the vehicle's coordinate system
           * the points in the simulator are connected by a Yellow line
           */

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          cout << msg <<endl;

          /* Latency
           * The purpose is to mimic real driving conditions where
           * the car does actuate the commands instantly.
           */
          long int latency_ms = (long int)(latency*1000);
          this_thread::sleep_for(chrono::milliseconds(latency_ms));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        /* Manual driving */
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  /* We don't need this since we're not using HTTP but if it's removed the
   * program doesn't compile :-(
   */
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t)
  {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req)
  {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length)
  {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
