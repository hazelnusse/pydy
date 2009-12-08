/* bicycle.h
 *
 * abstract data types Bicycle and Bicycle_benchmark
 *
 * Type Bicycle has 21 components, all parameters which describe the bicycle
 *
 * Type Bicycle_benchmark has 27 components, the parameters presented in
 * Meijaard et al., 2007.
 *
 */

typedef struct {
  double rr, rrt, rf, rft, lr, lf, ls;  // Wheel and frame geometric parameters
  double mr, mf;                  // Rear and front assembly mass
  double ICyy, IFyy;              // Out of plane wheel moments of inertia
  double IDxx, IDyy, IDzz, IDxz;  // Rear assembly inertia scalars
  double IExx, IEyy, IEzz, IExz;  // Front assembly inertia scalars
  double lrx, lrz, lfx, lfz;      // COM locations relative to wheel centers
} Bicycle_t;

typedef struct {
  double rr, rrt, rf, rft, w, c, lambda;  // Wheel and frame geometric parameters
  double mr, mb, mh, mf;          // R. whl, rear frame, fork, F. whl mass
  double IRxx, IRyy;              // Rear wheel inertia scalars, IRxx == IRzz
  double IBxx, IByy, IBzz, IBxz;  // Rear frame and rider inertia scalars
  double IHxx, IHyy, IHzz, IHxz;  // Front fork and handlebar inertia scalars
  double IFxx, IFyy;              // Front wheel inertia scalars, IFxx == IFzz
  double xb, zb, xh, zh;          // COM locations relative to rear contact
} Bicycle_benchmark_t;

/*
 * Converts Bicycle_t to Bicycle_benchmark_t
 */
void Bicycle_to_Bicycle_benchmark(Bicycle_benchmark_t *bike_b, const Bicycle_t *bike);

/*
 * Converts Bicycle_benchmark_t to Bicycle_t
 */
void Bicycle_benchmark_to_Bicycle(Bicycle_t *bike, const Bicycle_benchmark_t * bike_b);

/*
 * Calculates the time derivative of the eight generalized coordinates, given
 * the generalized coordinates, generalized speeds, and the bicycle parameters.
 */
void eval_qdots(double qd[], const double q[], const double u[], const Bicycle_t *bike);

/*
 * Calculates three dependent generalized speeds, given the coordinates (lean,
 * pitch, steer), the three independent generalized speeds, and the bicycle
 * parameters.  The indices of the independent speeds are specified by three
 * integer in the i arrays.  These integers must be distinct, and must take the
 * values of 0, 1, 2, 3, 4, or 5.
 */
void eval_dependent_speeds(double u[], double B[3][6], const int ind[], const int dep[]);


void mm3(double AB[3][3], double A[][3], double B[][3]);
void mv3(double Ab[3], double A[3][3], double b[3]);
void mdet3(double *det, double m[3][3]);
void madj3(double adj[3][3], double m[3][3]);
void form_constraint_matrix(double B[3][6], const double q[], const Bicycle_t *bike);


