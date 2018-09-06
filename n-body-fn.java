class nbody {
  final static double pi = 3.141592653589793;
  final static double solar_mass = 4.0 * pi * pi;
  final static double days_per_year = 365.24;

  static class Vector {
    double x, y, z;

    public Vector(double x, double y, double z) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
  }


  static Vector newVelocity(int bodyIdx, Vector velocity, Vector[] positions, double[] masses, double dt) {
    Vector position = positions[bodyIdx];
    double x = position.x;
    double y = position.y;
    double z = position.z;
    double dvx = 0.0;
    double dvy = 0.0;
    double dvz = 0.0;
    int len = positions.length;
    for (int i=0 ; i < len ; i++) {
      if (i != bodyIdx) {
        Vector p = positions[i];
        double dx = x - p.x;
        double dy = y - p.y;
        double dz = z - p.z;
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        double a_quantity = (masses[i] * dt) / (distance * distance * distance);
        dvx = dvx - dx * a_quantity;
        dvy = dvy - dy * a_quantity;
        dvz = dvz - dz * a_quantity;
      }
    }
    return new Vector(velocity.x + dvx, velocity.y + dvy, velocity.z + dvz);
  }


  static Vector[][] advance(Vector[] positions, Vector[] velocities, double[] masses, double dt) {
    int count = velocities.length;
    
    Vector[] newVelocities = new Vector[count];
    for (int i=0 ; i < count ; i++)
      newVelocities[i] = newVelocity(i, velocities[i], positions, masses, dt);

    Vector[] newPositions = new Vector[count];
    for (int i=0 ; i < count ; i++) {
      Vector position = positions[i];
      Vector velocity = newVelocities[i];
      newPositions[i] = new Vector(
        position.x + dt * velocity.x,
        position.y + dt * velocity.y,
        position.z + dt * velocity.z
      );
    }

    return new Vector[][] {newPositions, newVelocities};
  }


  static double energy(Vector[] positions, Vector[] velocities, double[] masses) {
    int count = positions.length;
    double e = 0.0;
    for (int i1=0 ; i1 < count ; i1++) {
      Vector p1 = positions[i1];
      Vector v1 = velocities[i1];
      double m1 = masses[i1];
      e = e + 0.5 * m1 * (v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
      for (int i2=i1+1 ; i2 < count ; i2++) {
        Vector p2 = positions[i2];
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        double dz = p1.z - p2.z;
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        e = e - (masses[i1] * masses[i2]) / distance;
      }
    }
    return e;
  }


  static void offsetMomentum(Vector[] velocities, double[] masses) {
    int count = velocities.length;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    for (int i=0 ; i < count ; i++) {
      Vector v = velocities[i];
      double m = masses[i];
      px = px + v.x * m;
      py = py + v.y * m;
      pz = pz + v.z * m;
    }
    Vector v0 = velocities[0];
    velocities[0] = new Vector(
      v0.x - px / solar_mass,
      v0.y - py / solar_mass,
      v0.z - pz / solar_mass
    );
  }


  static double[][] bodies = {
    // Sun
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, solar_mass},
    // Jupiter
    {  4.84143144246472090e00,
      -1.16032004402742839e00,
      -1.03622044471123109e-01,
       1.66007664274403694e-03 * days_per_year,
       7.69901118419740425e-03 * days_per_year,
      -6.90460016972063023e-05 * days_per_year,
       9.54791938424326609e-04 * solar_mass
    },
    // Saturn
    {  8.34336671824457987e00,
       4.12479856412430479e00,
      -4.03523417114321381e-01,
      -2.76742510726862411e-03 * days_per_year,
       4.99852801234917238e-03 * days_per_year,
       2.30417297573763929e-05 * days_per_year,
       2.85885980666130812e-04 * solar_mass
    },
    // Uranus
    {  1.28943695621391310e01,
      -1.51111514016986312e01,
      -2.23307578892655734e-01,
       2.96460137564761618e-03 * days_per_year,
       2.37847173959480950e-03 * days_per_year,
      -2.96589568540237556e-05 * days_per_year,
       4.36624404335156298e-05 * solar_mass
    },
    // Neptune
    {  1.53796971148509165e01,
      -2.59193146099879641e01,
       1.79258772950371181e-01,
       2.68067772490389322e-03 * days_per_year,
       1.62824170038242295e-03 * days_per_year,
      -9.51592254519715870e-05 * days_per_year,
       5.15138902046611451e-05 * solar_mass
    }
  };


  public static void main(String[] args) {
    int iteractions = Integer.parseInt(args[0]);

    int count = bodies.length;

    Vector[] positions = new Vector[count];
    Vector[] velocities = new Vector[count];
    double[] masses = new double[count];

    for (int i=0 ; i < count ; i++) {
      double[] body = bodies[i];
      positions[i] = new Vector(body[0], body[1], body[2]);
      velocities[i] = new Vector(body[3], body[4], body[5]);
      masses[i] = body[6];
    }

    offsetMomentum(velocities, masses);

    double e = energy(positions, velocities, masses);
    System.out.println(Double.toString(e));

    for (int i=0 ; i < iteractions ; i++) {
      Vector[][] posVel = advance(positions, velocities, masses, 0.01);
      positions = posVel[0];
      velocities = posVel[1];
    }

    e = energy(positions, velocities, masses);
    System.out.println(Double.toString(e));
  }
}