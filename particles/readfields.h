#pragma once

#include "vlsvreader2.h"
#include "vlsv_reader.h"
#include "vlsvreaderinterface.h"
#include "field.h"
#include <algorithm>
#include <vector>
#include <string>
#include <set>

template <class Reader>
class vlsvFieldReader {
protected:
	Reader r;
	std::vector<uint64_t> cellIds;
	Field E[2], B[2];
	std::string B_field_name;
	std::string E_field_name;
	std::string filename_pattern;
	int input_file_counter;

	/* Detect if volume-averaged fields are available */
	void detect_field_names() {

#ifdef DEBUG
		std::cerr << "Checking for volume-averaged fields... ";
#endif
		std::list<std::string> variableNames;
		std::string gridname("SpatialGrid");

		r.getVariableNames(gridname,variableNames);
		if(find(variableNames.begin(), variableNames.end(), std::string("B_vol"))!=variableNames.end()) {
#ifdef DEBUG
			std::cerr << "yep!" << std::endl;
#endif
			B_field_name = "B_vol";
			E_field_name = "E_vol";
		} else if(find(variableNames.begin(), variableNames.end(), std::string("B"))!=variableNames.end()) {
#ifdef DEBUG
			std::cerr << "Nope!" << std::endl;
#endif
			B_field_name = "B";
			E_field_name = "E";
		} else {
			std::cerr << "No B- or E-fields found! Strange file format?" << std::endl;
			exit(1);
		}
	}

public:

	/** Open filenames matching a given pattern
	 *  t0: Starting time
	 *  input_file_step: spacing between input steps (in seconds) */
	virtual void open(const std::string& _filename_pattern, double t0=0, double input_file_step=1.)=0;

	/* Read certain kinds of parameters */
	double readDoubleParameter(Reader& r, const char* name);
	/* Read a single-valued integer parameter */
	uint32_t readUintParameter(Reader& r, const char* name);

	/* Read the cellIDs into an array */
	std::vector<uint64_t> readCellIds(oldVlsv::Reader& r);
	std::vector<uint64_t> readCellIds(newVlsv::Reader& r);

	/* Read the "raw" field data in file order */
	std::vector<double> readFieldData(std::string& name, unsigned int numcomponents) {

			uint64_t arraySize=0;
			uint64_t vectorSize=0;
			uint64_t byteSize=0;
			vlsv::datatype::type dataType;
			std::list<std::pair<std::string,std::string> > attribs;
			attribs.push_back(std::pair<std::string,std::string>("name",name));
			if( r.getArrayInfo("VARIABLE",attribs, arraySize,vectorSize,dataType,byteSize) == false ) {
				std::cerr << "getArrayInfo returned false when trying to read VARIABLE \""
					<< name << "\"." << std::endl;
				exit(1);
			}

			if(dataType != vlsv::datatype::type::FLOAT || byteSize != 8 || vectorSize != numcomponents) {
				std::cerr << "Datatype of VARIABLE \"" << name << "\" entries is not double." << std::endl;
				exit(1);
			}

			/* Allocate memory for the data */
			std::vector<double> buffer(arraySize*vectorSize);

			if( r.readArray("VARIABLE",attribs,0,arraySize,(char*) buffer.data()) == false) {
				std::cerr << "readArray faied when trying to read VARIABLE \"" << name << "\"." << std::endl;
				exit(1);
			}

			return buffer;
		}

	/* Read a field (this relies on the cellIDs having been read already) */
	Field readField(std::string name, unsigned int numcomponents) {
		Field f;
		double min[3], max[3], time;
		uint64_t cells[3];
		min[0] = readDoubleParameter(r,"xmin");
		min[1] = readDoubleParameter(r,"ymin");
		min[2] = readDoubleParameter(r,"zmin");
		max[0] = readDoubleParameter(r,"xmax");
		max[1] = readDoubleParameter(r,"ymax");
		max[2] = readDoubleParameter(r,"zmax");
		cells[0] = readUintParameter(r,"xcells_ini");
		cells[1] = readUintParameter(r,"ycells_ini");
		cells[2] = readUintParameter(r,"zcells_ini");
		time = readDoubleParameter(r,"t");

		f.data.resize(cells[0]*cells[1]*cells[2]*4);

		std::vector<double> buffer = readFieldData(name,numcomponents);

		/* Sanity-check stored data sizes */
		if(numcomponents*cellIds.size() != buffer.size()) {
			std::cerr << numcomponents << " * cellIDs.size (" << cellIds.size() 
				  << ") != buffer.size (" << buffer.size() << ") while reading "
				  << name << "!" << std::endl;
			exit(1);
		}

		/* Set field size */
		for(int i=0; i<3;i++) {
			/* Volume-centered values -> shift by half a cell in all directions*/
			f.dx[i] = (max[i]-min[i])/cells[i];
			double shift = f.dx[i]/2;
			f.min[i] = min[i]+shift;
			f.max[i] =  max[i]+shift;
			f.cells[i] = cells[i];
		}
		f.time = time;

		/* So, now we've got the cellIDs, the mesh size and the field values,
		 * we can sort them into place */
		for(uint i=0; i< cellIds.size(); i++) {
			uint64_t c = cellIds[i];
			int64_t x = c % cells[0];
			int64_t y = (c /cells[0]) % cells[1];
			int64_t z = c /(cells[0]*cells[1]);

			double* ftgt = f.getCellRef(x,y,z);
			ftgt[0] = buffer[3*i];
			ftgt[1] = buffer[3*i+1];
			ftgt[2] = buffer[3*i+2];
		}

		return f;
	}

	Field read_velocity_field() {

		Field V = readField("rho_v",3u);
		Field rho = readField("rho",1u);

		for(uint i=0; i< cellIds.size(); i++) {
			V.data[4*i] /= rho.data[4*i];
			V.data[4*i+1] /= rho.data[4*i];
			V.data[4*i+2] /= rho.data[4*i];
		}

		return V;
	}

	Interpolated_Field interpolate_B(double time) {
		return Interpolated_Field(B[0],B[1],time);
	}
	Interpolated_Field interpolate_E(double time) {
		return Interpolated_Field(B[0],B[1],time);
	}

	/** Open the file data required for the next timestep.
	 *  If the corresponding file is already opened, simply do nothing.
	 *  Return value specifies whether a new file was opened. */
	virtual bool open_next_timestep(double t) = 0;
	virtual bool read_next_timestep(double t, Field& E0, Field& E1,
			Field& B0, Field& B1) {

		if(open_next_timestep(t)) {
			uint64_t cells[3];
			cells[0] = readUintParameter(r,"xcells_ini");
			cells[1] = readUintParameter(r,"ycells_ini");
			cells[2] = readUintParameter(r,"zcells_ini");

			/* Read CellIDs and Field data */
			cellIds = readCellIds(r);
			std::string name(B_field_name);
			std::vector<double> Bbuffer = readFieldData(name,3u);
			name = E_field_name;
			std::vector<double> Ebuffer = readFieldData(name,3u);

			/* Assign them, without sanity checking */
			/* TODO: Is this actually a good idea? */
			for(uint i=0; i< cellIds.size(); i++) {
				uint64_t c = cellIds[i];
				int64_t x = c % cells[0];
				int64_t y = (c /cells[0]) % cells[1];
				int64_t z = c /(cells[0]*cells[1]);

				double* Etgt = E1.getCellRef(x,y,z);
				double* Btgt = B1.getCellRef(x,y,z);
				Etgt[0] = Ebuffer[3*i];
				Etgt[1] = Ebuffer[3*i+1];
				Etgt[2] = Ebuffer[3*i+2];
				Btgt[0] = Bbuffer[3*i];
				Btgt[1] = Bbuffer[3*i+1];
				Btgt[2] = Bbuffer[3*i+2];
			}

			return true;
		}

		return false;
	}

};

template <class Reader>
class vlsvForwardFieldReader: public vlsvFieldReader<Reader> {
protected:
	// Darn you, C++ template inheritance
	using vlsvFieldReader<Reader>::E;
	using vlsvFieldReader<Reader>::B;
	using vlsvFieldReader<Reader>::cellIds;
	using vlsvFieldReader<Reader>::E_field_name;
	using vlsvFieldReader<Reader>::B_field_name;
	using vlsvFieldReader<Reader>::r;
	using vlsvFieldReader<Reader>::input_file_counter;
	using vlsvFieldReader<Reader>::filename_pattern;
	using vlsvFieldReader<Reader>::readFieldData;
	using vlsvFieldReader<Reader>::detect_field_names;

public:
	virtual void open(const std::string& _filename_pattern, double t0=0, double input_file_step=1.) {
		char filename_buffer[256];
		filename_pattern = _filename_pattern;
		input_file_counter = floor(t0/input_file_step)-1;
		snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);


#ifdef DEBUG
		std::cerr << "Opening " << filename_buffer << "...";
#endif
		r.open(filename_buffer);
#ifdef DEBUG
		std::cerr <<"ok." << std::endl;
#endif
		detect_field_names();

		/* Set dummy times so that the next timestep will always be read */
		E[0].time = 99999999999999.;
		E[1].time = -99999999999999.;

		/* Read the cellIDs so they can next be used to sort field data
                 * into place */
		/* Read the MESH, yielding the CellIDs */
		cellIds = readCellIds(r);

	}
	virtual bool open_next_timestep(double t) {
		char filename_buffer[256];
		bool retval = false;

		while(t < E[0].time || t>= E[1].time) {
#ifdef DEBUG
			std::cerr << "Reading next timestep. t = " << t << "  E[0] is at " <<
				E[0].time << ",  E[1] is at " << E[1].time << std::endl;
#endif
			snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);
			input_file_counter += 1;

			r.close();
#ifdef DEBUG
		std::cerr << "Opening " << filename_buffer << "...";
#endif
			r.open(filename_buffer);
#ifdef DEBUG
		std::cerr <<"ok." << std::endl;
#endif
			retval = true;
			double new_t = readDoubleParameter(r,"t");
			cellIds = readCellIds(r);

			if(new_t <= t) {
				/* If the other timestep is already initialized, push onwards */
				if(E[0].time < 99999999999999.) {
					E[1] = E[0];
					B[1] = B[0];
				}
				E[0] = readField(E_field_name,3);
				B[0] = readField(B_field_name,3);
			}
			if(new_t > t) {
				if(E[1].time > -99999999999999.) {
					E[0] = E[1];
					B[0] = B[1];
				}
				E[1] = readField(E_field_name,3);
				B[1] = readField(B_field_name,3);
			}
		}
		return retval;
	}



};

template <class Reader>
class vlsvBackwardFieldReader: public vlsvFieldReader<Reader> {
protected:
	// Darn you, C++ template inheritance
	using vlsvFieldReader<Reader>::E;
	using vlsvFieldReader<Reader>::B;
	using vlsvFieldReader<Reader>::cellIds;
	using vlsvFieldReader<Reader>::E_field_name;
	using vlsvFieldReader<Reader>::B_field_name;
	using vlsvFieldReader<Reader>::r;
	using vlsvFieldReader<Reader>::input_file_counter;
	using vlsvFieldReader<Reader>::filename_pattern;
	using vlsvFieldReader<Reader>::readFieldData;
	using vlsvFieldReader<Reader>::detect_field_names;
public:
	virtual void open(const std::string& _filename_pattern, double t0=0, double input_file_step=1.) {
		char filename_buffer[256];
		filename_pattern = _filename_pattern;
		input_file_counter = floor(t0/input_file_step)+1;
		snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);


#ifdef DEBUG
		std::cerr << "Opening " << filename_buffer << "...";
#endif
		r.open(filename_buffer);
#ifdef DEBUG
		std::cerr <<"ok." << std::endl;
#endif
		detect_field_names();

		/* Set dummy times so that the next timestep will always be read */
		E[0].time = 99999999999999.;
		E[1].time = -99999999999999.;

		/* Read the cellIDs so they can next be used to sort field data
                 * into place */
		/* Read the MESH, yielding the CellIDs */
		cellIds = readCellIds(r);

	}
	virtual bool open_next_timestep(double t) {
		char filename_buffer[256];
		bool retval = false;

		while(t < E[0].time || t>= E[1].time) {
#ifdef DEBUG
			std::cerr << "Read next timestep. t = " << t << "  E[0] is at " <<
				E[0].time << ",  E[1] is at " << E[1].time << std::endl;
#endif /* DEBUG */
			snprintf(filename_buffer,256,filename_pattern.c_str(),input_file_counter);
			input_file_counter -= 1;

			r.close();
#ifdef DEBUG
		std::cerr << "Opening " << filename_buffer << "...";
#endif
			r.open(filename_buffer);
#ifdef DEBUG
		std::cerr <<"ok." << std::endl;
#endif
			retval = true;
			double new_t = readDoubleParameter(r,"t");
			cellIds = readCellIds(r);

			if(new_t <= t) {
				/* If the other timestep is already initialized, push onwards */
				if(E[0].time < 99999999999999.) {
					E[1] = E[0];
					B[1] = B[0];
				}
				E[0] = readField(E_field_name,3);
				B[0] = readField(B_field_name,3);
			}
			if(new_t > t) {
				if(E[1].time > -99999999999999.) {
					E[0] = E[1];
					B[0] = B[1];
				}
				E[1] = readField(E_field_name,3);
				B[1] = readField(B_field_name,3);
			}
		}
		return retval;
	}
};


/* For debugging purposes - dump a field into a png file */
void debug_output(Field& F, const char* filename);
