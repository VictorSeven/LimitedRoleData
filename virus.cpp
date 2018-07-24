#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

//Include network code
#include "CNet.cpp"

//Define global quantities
#define L 100
#define N L*L

#define SUS 0
#define INF 1
#define REC 2
#define EMP 3

using namespace std;

//Declaration of functions (explained on implementation below)
void initial_conditions(CNetwork<int> &lattice, double p, double p_e, double beta0, double ti0, mt19937 &gen, uniform_real_distribution<double> &ran_u);
void read_initial_conditions(CNetwork<int> &lattice, int savefile, string file_name);
inline void write_videos(int iteration, int &video, int its, int frames, char* file_name, CNetwork<int> &lattice, ofstream &output_v, ofstream &output_b, ofstream &output_t, ofstream &output_time);
void update_lattice(int iteration, CNetwork<int> &lattice, CNetwork<int> &old_lattice, double tr, double mu, double dt, double dbeta, double dti, double elapsed_time, int &total_infecciones,
                    mt19937& gen, uniform_real_distribution<double> &ran_u, ofstream &output);
void update_part_lattice(int start, int end, CNetwork<int> &lattice, CNetwork<int> &old_lattice, double tr, double mu, double dt, double dbeta, double dti, double elapsed_time,
							double &mean_beta, double &mean_beta2, double &mean_ti, double &mean_ti2, int &n_sus, int &n_inf, int &n_rec, mt19937& gen, uniform_real_distribution<double> &ran_u);
void rotate_particles(CNetwork<int> &lattice, bool &start_0, int D, mt19937& gen, uniform_real_distribution<double> &ran_u);
inline void overwrite_values(CNetwork<int> &lattice, int old_n, int new_n);
vector<double> correlation_length(CNetwork<int> &lattice);
int dist(int i1, int j1, int i2, int j2);
inline int isingmap(int type);

// ====================================================================================== //
// ************************************************************************************** //
// ***************************           MAIN CODE          ***************************** //
// ************************************************************************************** //
// ====================================================================================== //



int main(int argc, char* argv[])
{
    int i,j,k; //Counters
    const double tr = 1.0; //Tau_r -> constant time scale

    //Used to take out data
    ofstream output;
    ofstream output_v, output_b, output_t, output_time;
    ofstream output_corr;

    int random_seed; //Seed of the RNG
    char* file_name; //Where to save file

    int recover;

    vector<double> corr_len; //Stores correlation length

    int its;        //Number of iterations
    int frames;     //Number of frames (for video outputs)

    int video;

    double mu = 0.01;  //Mutation rate

    //Evolutionary parameters
    double beta = 0.6;
    double ti = 1.5;

    //Time
    double elapsed_time = 0.0;

    //Change of quantities (values as indicated by Ballegooijen&Boerlijst)
    double dt = 0.01;
    double dbeta = 0.01;
    double dti = 0.01;

	bool start_0 = true; //True = starts rotation in (0,0); False = starts rotation in (1,1)
	int D=1; //Number of times we mix (Fixed to 1 in ALL experiments)
    int D_counter = 10;     //This avoid doing a mix in every iteration. Controls the diffusion constant d=D/D_counter=1/D_counter

    //Counter of total infections
    int total_infecciones = 0;

    //Synchronous update needs two lattices
    CNetwork<int> lattice(N, false); //L^2 nodes, unweighted network that store a int in each node
    CNetwork<int> old_lattice(N, false);

    //Make square lattice with 8 neighbours
    lattice.create_lattice_8(L);
    //Define properties for each node
    lattice.define_property("elapsed_time", "double", true);
    lattice.define_property("beta", "double", true);
    lattice.define_property("ti", "double", true);



    //Read arguments from console
    if (argc < 9)
    {
        cout << "Error: must include its, beta, ti, random, D and file_name parameters" << endl;
        return 0;
    }
    else
    {
        its = stoi(argv[1]);
        beta = stod(argv[2]);
        ti = stod(argv[3]);
        random_seed = stoi(argv[4]);
        D_counter = stoi(argv[5]);
        file_name = argv[6];
        frames = stoi(argv[7]);
        recover = stoi(argv[8]);
    }

    //Init the RNG
    mt19937 gen(random_seed);
    uniform_real_distribution<double> ran_u(0.0,1.0);


    elapsed_time = 0.0; //Init elapsed time

    //Allow user to recover an existent simulation from any point
    if (recover <= 0)
    {
        initial_conditions(lattice, 0.05, 0.0, beta, ti, gen, ran_u);
        video = 0;
    }
    else
    {
        read_initial_conditions(lattice, recover, string(file_name));
        video = recover+1;
    }


    //Init this variable
    old_lattice = lattice;


    //Save results here
    output.open("output/"+string(file_name));
    output_corr.open("output/corr_"+string(file_name));

    //Main MonteCarlo (MC) loop, where dynamics run
    for (i=0; i < its; i++)
    {
        //Do another MC step
        update_lattice(i, lattice, old_lattice, tr, mu, dt, dbeta, dti, elapsed_time, total_infecciones, gen, ran_u, output);

        //Shuffle the system
        //UNCOMMENT THIS FOR NON-ZERO DIFFUSION
        //if (i % D_counter == 0) rotate_particles(lattice, start_0, D, gen, ran_u);

        //Write a fun video
        //UNCOMMENT THIS TO SEE VIDEO OF DYNAMICS
        //write_videos(i, video, its, frames, file_name, lattice, output_v, output_b, output_t, output_time);

        //Compute the correlation lenght
        //UNCOMMENT THIS TO COMPUTE CORRELATION LENGHT
        /*
        if (i % int(its/frames) == 0)
        {
            corr_len = correlation_length(lattice);
            for (j=0; j < corr_len.size(); j++)
            {
                output_corr << corr_len[j] << " ";
            }
            output_corr << endl;
        }
        */

        //Update the lattice and increase time
        old_lattice = lattice;

        elapsed_time += dt;
    }

    //Close the outputs
    output_corr.close();
    output.close();

    return 0;
}




// ====================================================================================== //
// ************************************************************************************** //
// ***************************    INITIAL CONDITIONS        ***************************** //
// ************************************************************************************** //
// ====================================================================================== //

/** This function sets the initial conditions for the lattice. This must be invoked before starting the simulation.
* \param lattice The CNetwork that represents the system
* \param p Probability of having an infected site
* \param p_e Probability of having an empty site
* \param beta0 Initial value of beta
* \param ti0 Initial value of ti0
* \param gen RNG
* \param ran_u Random numbers in the interval [0,1]
*/
void initial_conditions(CNetwork<int> &lattice, double p, double p_e, double beta0, double ti0, mt19937 &gen, uniform_real_distribution<double> &ran_u)
{
    int i;
    double random;
    //For each cell, set the values
    for (i=0; i < N; i++)
    {
        random = ran_u(gen);
        if (random <= p) lattice.set_value(i, INF); //Infected
        else if (random <= p + p_e) lattice.set_value(i, EMP); //Empty cells
        else lattice.set_value(i, SUS);

        lattice.set_value("beta", i, beta0);
        lattice.set_value("ti", i, ti0);
        lattice.set_value("elapsed_time", i, 0.0);
    }


}

/** This function is able to recover a simulation that started before and was saved. This is useful to do
* long simulations without worrying of any data loss in the cluster, or to decide if it worth to continue a computation or not.
* In order to use this function the writing of videos must be UNCOMMENTed.
* \param lattice The CNetwork that represents the system
* \param savefile Index of the savefile, which is reprsented by the "frame" variable when writing videos.
* \param file_name The name of the file
*/
void read_initial_conditions(CNetwork<int> &lattice, int savefile, string file_name)
{
    ifstream input_state, input_beta, input_ti, input_time;

    int i,j;

    int infected;
    double beta, ti, time;

    //Open the four files that store the matrix
    input_state.open("videos/video_" + to_string(savefile) + string(file_name));
    input_beta.open("videos/video_beta_" + to_string(savefile) + string(file_name));
    input_ti.open("videos/video_ti_" + to_string(savefile) + string(file_name));
    input_time.open("videos/video_time_" + to_string(savefile) + string(file_name));

    //Read info and load it into the variable
    for (i=0; i < L; i++)
    {
        for (j=0; j < L; j++)
        {
            //Read data from each file
            input_state >> infected;
            input_beta >> beta;
            input_ti >> ti;
            input_time >> time;

            //Take it into lattice
            lattice.set_value(i+j*L, infected);
            lattice.set_value("beta", i+j*L, beta);
            lattice.set_value("ti", i+j*L, ti);
            lattice.set_value("elapsed_time", i+j*L, time);
        }
    }
    input_state.close();
    input_beta.close();
    input_ti.close();
    input_time.close();

    return;
}


/** This function writes a new frame of a video. It is also used to make checkpoints, in order to store
* the current state of the simulation, which may be recovered after.
* \param iteration Current iteration of main MC loop
* \param video Frame index
* \param its Total number of iterations
* \param frames Total number of frames we desire
* \param file_name The name of the file
* \param lattice The CNetwork that represents the lattice
* \param output_v Represents the file where the state of nodes will be saved
* \param output_b Represents the file where the value of each beta_i will be saved
* \param output_t Represents the file where the value of each tauI_i will be saved
* \param output_time Represents the file where  internal time counter of nodes will be saved
*/
void write_videos(int iteration, int &video, int its, int frames, char* file_name, CNetwork<int> &lattice, ofstream &output_v, ofstream &output_b, ofstream &output_t, ofstream &output_time)
{
    int j,k;

    //In order to have the specified frames, save a frame each its/frames steps
    if (iteration % int(its / frames) == 0)
    {
        //Open the corresponding folders
        output_v.open("videos/video_" + to_string(video) + string(file_name));
        output_b.open("videos/video_beta_" + to_string(video) + string(file_name));
        output_t.open("videos/video_ti_" + to_string(video) + string(file_name));
        output_time.open("videos/video_time_" + to_string(video) + string(file_name));

        //Save all states
        for (j=0; j < L; j++)
        {
            for (k=0; k < L; k++)
            {
                output_v << lattice.get_value(k+j*L) << " "; //Write who is infected
                output_time << lattice.get_value_d("elapsed_time", k+j*L); //Store internal time to re-do process again

                //For infected, also store information about the strain:
                if (lattice.get_value(k+j*L) == INF)
                {
                    output_b << lattice.get_value_d("beta", k+j*L) << " ";
                    output_t << lattice.get_value_d("ti", k+j*L) << " ";
                }
                else
                {
                    output_b << -1.0 << " ";
                    output_t << -1.0 << " ";
                }

            }
            output_v << endl;
            output_b << endl;
            output_t << endl;
            output_time << endl;
        }

        output_v.close();
        output_b.close();
        output_t.close();
        output_time.close();

        video++; //Increase the counter
    }
}



// ====================================================================================== //
// ************************************************************************************** //
// ***************************    ALGORITHM IMPLEMENTATION  ***************************** //
// ************************************************************************************** //
// ====================================================================================== //




/** This function implements the algorithm discussed in the paper, and it is the main responsible for the dynamics of the system.
* \param iteration Current iteration of main MC loop
* \param lattice The CNetwork that represents the lattice. This will be modified.
* \param old_lattice The CNetwork in the last time step
* \param tr Tau_r. Fixed to 1 in most experiments
* \param mu Mutation rate
* \param dt Timestep for the algorithm
* \param dbeta Change in beta each mutation
* \param dti Change in tauI each mutation
* \param elapsed_time Time of the simulated experiment
* \param total_infecciones Total number of infections
* \param gen RNG
* \param ran_u Random number between 0 and 1
* \param output Represents the file where the results of averages will be saved
*/
void update_lattice(int iteration, CNetwork<int> &lattice, CNetwork<int> &old_lattice, double tr, double mu, double dt, double dbeta, double dti, double elapsed_time, int &total_infecciones,
                    mt19937& gen, uniform_real_distribution<double> &ran_u, ofstream &output)
{
    int i,j; //Counters

    int cell, neigh; //Indices of selected cell and its neighbours

    //Values of these cells
    double beta_neigh;
    double beta_cell;
    double ti_cell;

    bool finish;

    //Auxiliary variable
    double aux;

    //Total number of each kind
    int n_sus, n_inf, n_rec;

    //In order to compute averages for first and second moments
    double mean_beta, mean_ti;
    double mean_beta2, mean_ti2;
    double varbeta, varti;


    double p_infection; //Probability of getting infected
    double beta_infection; //Sum of beta of neigh
    vector<double> p_neigh_inf = vector<double>(0); //Probabiliy of infection of each neighbour
    double p_neigh_inf_total; //Sum of probabilities of neighbours
    double right; //To select infection neighbour
    int n_infected; //Number of infected neighbours
    vector<int> who_infect = vector<int>(0); //Who to infect

    //Init everything
    n_sus = n_inf = n_rec = 0;
    mean_beta = mean_beta2 = mean_ti = mean_ti2 = 0.0;

    //Update all the lattice
    for (i = 0; i < N; i++)
    {
        //Get the current cell and its values
        cell = old_lattice.get_value(i);
        beta_cell = old_lattice.get_value_d("beta", i);
        ti_cell = old_lattice.get_value_d("ti", i);

        //Update elapsed time
        lattice.set_value("elapsed_time", i, old_lattice.get_value_d("elapsed_time", i) + dt);


        if (cell == SUS)
        {
            n_sus++;

            //Prepare variables in order see if I should become infected
            j = 0;
            beta_infection = 0.0;
            n_infected = 0;
            p_neigh_inf_total = 0.0;
            who_infect.clear();
            p_neigh_inf.clear();

            //*************************************  IMPORTANT  *******************************************//
            //IThe reviewer reminded me that when Poisson process compete, we can check who was
            //the one infecting with a simple probability. This modification was included in the algorithm
            //explained in the paper (it is the Taylor expansion in step 2a), but the code was not modified.
            //Including the modification may speed up the code in a significant amount!
            //*********************************************************************************************//

            //Loop over neighs
            while (j < old_lattice.get_num_neighs(i))
            {
                //Get neigh
                neigh = old_lattice.get_neigh_at(i, j);

                //If it is infected, check if the infection transmits
                if (old_lattice.get_value(neigh) == INF)
                {

                    beta_neigh = old_lattice.get_value_d("beta", neigh); //Infectiveness of the neighbour
                    beta_infection += beta_neigh;


                    p_neigh_inf.push_back(1.0 - exp(- beta_neigh * dt)); //Store probability of infection
                    p_neigh_inf_total += p_neigh_inf[n_infected]; //Sum it, will be useful later
                    who_infect.push_back(neigh);

                    n_infected++;
                }
                j++;
            } //End while


            //Probability of being infected. Draw a random number
            p_infection = 1.0 - exp(- beta_infection * dt);
            aux = ran_u(gen);

            //If we get infected, then select who is the infected
            if (aux <= p_infection)
            {
                aux /= p_infection; //Put aux between 0 and 1 to select again in the rest of the line

                //Select on the line who is going to be the infector:
                neigh = 0;
                right = p_neigh_inf[0] / p_neigh_inf_total;

                while (aux > right)
                {
                    neigh++;
                    right += p_neigh_inf[neigh] / p_neigh_inf_total; //Get probability of infection of next neighbour, rescaled from 0 to 1
                }
                neigh = who_infect[neigh]; //Replace the index with the one of the lattice

                lattice.set_value(i, INF); //Infect this node
                lattice.set_value("beta", i, old_lattice.get_value_d("beta", neigh)); //Get the new beta of the neighbour
                lattice.set_value("ti", i, old_lattice.get_value_d("ti", neigh)); //Get the new tau_i of the neighbour
                lattice.set_value("elapsed_time", i, 0.0); //Set the elapsed time

            }


        }
        else if (cell == INF) //and lattice.get_value_d("elapsed_time", i) > ti)
        {
            n_inf++;
            mean_beta += beta_cell;
            mean_beta2 += beta_cell * beta_cell;
            mean_ti += ti_cell;
            mean_ti2 += ti_cell * ti_cell;

            if (lattice.get_value_d("elapsed_time", i) >= ti_cell)
            {
                lattice.set_value(i, REC);
                lattice.set_value("elapsed_time", i, 0.0);
            }

            //Mutate infected individuals
            //UNCOMMENT TO ACTIVATE MUTATION
            aux = ran_u(gen);
            if (aux <= mu * dt)
            {
                aux /= mu * dt; //Put aux between 0 and 1
                if (aux <= 0.25)
                {
                    lattice.set_value("beta", i,  beta_cell + dbeta);
                    lattice.set_value("ti", i, ti_cell + dti);
                }
                else if (aux <= 0.5)
                {
                    lattice.set_value("beta", i,  beta_cell + dbeta);
                    lattice.set_value("ti", i, ti_cell - dti);
                }
                else if (aux <= 0.75)
                {
                    lattice.set_value("beta", i,  beta_cell - dbeta);
                    lattice.set_value("ti", i, ti_cell + dti);
                }
                else
                {
                    lattice.set_value("beta", i,  beta_cell - dbeta);
                    lattice.set_value("ti", i, ti_cell - dti);
                }

            }


        }
        else if (cell == REC)
        {
            n_rec++;
            if (lattice.get_value_d("elapsed_time", i) >= tr)
            {
                lattice.set_value(i, SUS);
                lattice.set_value("elapsed_time", i, 0.0);
            }
        }
        else
        {
            ; //Empty cells do nothing
        }

    }


    total_infecciones += n_inf; //Increase number of total infections between two writes


    //Finish means and compute variances
    if (n_inf != 0)
    {
        mean_beta /= 1.0 * n_inf;
        mean_ti /= 1.0 * n_inf;
        varbeta = mean_beta2 / (1.0 * n_inf) - mean_beta * mean_beta;
        varti = mean_ti2 / (1.0 * n_inf) - mean_ti * mean_ti;
    }
    else
    {
        mean_beta = 0.0;
        mean_ti = 0.0;
        varbeta =  -1.0;
        varti = -1.0;
    }


    //One of each 1000 steps, save a data point!
    if (iteration % 1000 == 0)
    {
        output << elapsed_time << " " << mean_beta << " " << mean_ti << " " << total_infecciones/(1.0*N) << endl;
        total_infecciones = 0;
    }

}



/** This function implements the Toffolis algorithm for diffusion.
* \param lattice The CNetwork that represents the lattice. This will be modified.
* \param start_0 If we start the mixing at (0,0) or at (1,1)
* \param D Number of times we mix
* \param gen RNG
* \param ran_u Random number between 0 and 1
*/
void rotate_particles(CNetwork<int> &lattice, bool &start_0, int D, mt19937& gen, uniform_real_distribution<double> &ran_u)
{
	int j; //Counter
	int x,y; //Coordinates of each cell

	int c, cr, cu, cur; //Cell, Cell-Right, Cell-Up, Cell-Up-Right

	//To copy all properties of a node
	int copia_tipo;
	double copia_beta, copia_ti, copia_internal;

    //Mix D times (D=1 throught all simulations)
	for (j=0; j < D; j++)
	{
	    //Start mixing from (0,0)
		if (start_0)
		{
			if (ran_u(gen) < 0.5) //Move everything clockwise
			{

				y = 0;
				while (y < L)
				{
					x = 0;
					while(x < L)
					{
						//Get coordinates of the cells in a 2x2 square
						c = x + y * L;
						cr = c + 1;
						cu = c + L;
						cur = cu + 1;

						//Copy the first one here
						copia_tipo = lattice.get_value(c);
						copia_beta = lattice.get_value_d("beta", c);
						copia_ti = lattice.get_value_d("ti", c);
						copia_internal = lattice.get_value_d("elapsed_time", c);

						overwrite_values(lattice, c, cr);  //Get values from the one to the right and copy here
						overwrite_values(lattice, cr, cur); //Now overwrite the ones at the right with the ones up
						overwrite_values(lattice, cur, cu); //Overwrite the right up corner with the up left...

						//And use the copy to recover this the up left
						lattice.set_value(cu, copia_tipo);
						lattice.set_value("beta", cu, copia_beta);
						lattice.set_value("ti", cu, copia_ti);
						lattice.set_value("elapsed_time", cu, copia_internal);
						x += 2;
					}
					y += 2;
				}
			}
			else //Move everything counter-clockwise
			{
				y = 0;
				while (y < L)
				{
					x = 0;
					while(x < L)
					{
						//Get coordinates of the cells in a 2x2 square
						c = x + y * L;
						cr = c + 1;
						cu = c + L;
						cur = cu + 1;

                        //Copy the ones here
						copia_tipo = lattice.get_value(c);
						copia_beta = lattice.get_value_d("beta", c);
						copia_ti = lattice.get_value_d("ti", c);
						copia_internal = lattice.get_value_d("elapsed_time", c);

						overwrite_values(lattice, c, cu);
						overwrite_values(lattice, cu, cur);
						overwrite_values(lattice, cur, cr);

						lattice.set_value(cr, copia_tipo);
						lattice.set_value("beta", cr, copia_beta);
						lattice.set_value("ti", cr,  copia_ti);
						lattice.set_value("elapsed_time", cr, copia_internal);

						x += 2;
					}
					y += 2;
				}

			}
		}
		else //Start mixing from (1,1)
		{
			if (ran_u(gen) < 0.5) //Move everything clockwise
			{

				y = 1;
				while (y < L)
				{
					x = 1;
					while(x < L)
					{
						//Get coordinates of the cells in a 2x2 square -remembering periodic conditions
                        c = x + y * L; //Maximum of c is x=L-1
                        cr = ((x + 1) % L) + y * L;
                        cu = x + ((y + 1) % L) * L;
                        cur = ((x + 1) % L) + ((y + 1) % L) * L;

						//Copy the first one here
						copia_tipo = lattice.get_value(c);
						copia_beta = lattice.get_value_d("beta", c);
						copia_ti = lattice.get_value_d("ti", c);
						copia_internal = lattice.get_value_d("elapsed_time", c);

						overwrite_values(lattice, c, cr);  //Get values from the one to the right and copy here
						overwrite_values(lattice, cr, cur); //Now overwrite the ones at the right with the ones up
						overwrite_values(lattice, cur, cu); //Overwrite the right up corner with the up left...

						//And use the copy to recover this the up left
						lattice.set_value(cu, copia_tipo);
						lattice.set_value("beta", cu, copia_beta);
						lattice.set_value("ti", cu, copia_ti);
						lattice.set_value("elapsed_time", cu, copia_internal);
						x += 2;
					}
					y += 2;
				}
			}
			else //Move everything counterclockwise
			{
				y = 1;
				while (y < L)
				{
					x = 1;
					while(x < L)
					{
						//Get coordinates of the cells in a 2x2 square
                        c = x + y * L; //Maximum of c is x=L-1
                        cr = ((x + 1) % L) + y * L;
                        cu = x + ((y + 1) % L) * L;
                        cur = ((x + 1) % L) + ((y + 1) % L) * L;

						copia_tipo = lattice.get_value(c);
						copia_beta = lattice.get_value_d("beta", c);
						copia_ti = lattice.get_value_d("ti", c);
						copia_internal = lattice.get_value_d("elapsed_time", c);

						overwrite_values(lattice, c, cu);
						overwrite_values(lattice, cu, cur);
						overwrite_values(lattice, cur, cr);

						lattice.set_value(cr, copia_tipo);
						lattice.set_value("beta", cr, copia_beta);
						lattice.set_value("ti", cr,  copia_ti);
						lattice.set_value("elapsed_time", cr, copia_internal);

						x += 2;
					}
					y += 2;
				}

			}
		}//End start_0 else

		start_0 = not start_0; //Next time we will have the change in the other direction.


	} //End for

}



// ====================================================================================== //
// ************************************************************************************** //
// ***************************      AUXILIARY FUNCTIONS     ***************************** //
// ************************************************************************************** //
// ====================================================================================== //


/** This function takes the properties of the node NEW_N and copy them into the node OLD_N.
* It is used in order to reduce the number of code in the rotate function.
* \param lattice The CNetwork that represents the lattice. This will be modified.
* \param old_n Node to overwrite
* \param new_n Source node to copy
*/
inline void overwrite_values(CNetwork<int> &lattice, int old_n, int new_n)
{
	lattice.set_value(old_n, lattice.get_value(new_n) );
	lattice.set_value("beta", old_n, lattice.get_value_d("beta", new_n) );
	lattice.set_value("ti", old_n, lattice.get_value_d("ti", new_n));
	lattice.set_value("elapsed_time", old_n, lattice.get_value_d("elapsed_time", new_n));
}


/** Computes in a simple way the correlantion lenght of the given lattice
* \param lattice The CNetwork that represents the lattice. This will be modified.
*/
vector<double> correlation_length(const CNetwork<int> &lattice)
{
    int i,j;
    int i1,j1,i2,j2;
    vector<double> corr = vector<double>(1+L/2, 0.0); //Due to the definition of distance, maximum distance is L/2...
    vector<int> ncount = vector<int>(1+L/2, 0); //To count how many elements I have summed
    int d;

    //For all lattice
    for (i=0; i < N; i++)
    {
        //Get (x,y) coordinates of node i
        j1 = int(floor(i / L));
        i1 = i - j1 * L;

        //Now, for the rest of the lattice...
        for (j=0; j < N; j++)
        {
            //Get (x,y) coordinates of j
            j2 = int(floor(j / L));
            i2 = j - j2 * L;

            //Compute the distance and check the correlation at this distance
            d = dist(i1, j1, i2, j2);
            corr[d] += isingmap(lattice.get_value(i)) * isingmap(lattice.get_value(j));
            ncount[d]++;
        }
    }

    //Finish normalizing by the number of times a distance is found
    for (i=0; i < corr.size(); i++)
    {
        corr[i] /= ncount[i];
    }

    return corr;
}

//Computes distance between two nodes in the 8-neighbours lattice
inline int dist(int i1, int j1, int i2, int j2)
{
    int dx = abs(i1 - i2);
    int dy = abs(j1 - j2);

    dx = dx > L/2 ? L - dx : dx;
    dy = dy > L/2 ? L - dy : dy;

    return max(dx, dy);
}

//+1 to infected, -1 to the others -so we can use the definition of correlation of the Ising model
inline int isingmap(int type)
{
    if (type == INF) return 1;
    else return -1;
}
