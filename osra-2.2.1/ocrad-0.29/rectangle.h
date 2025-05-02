/* GNU Ocrad - Optical Character Recognition program
   Copyright (C) 2003-2024 Antonio Diaz Diaz.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

class Rectangle
  {
  int left_, top_, right_, bottom_;		// inclusive coordinates

public:
  Rectangle( const int l, const int t, const int r, const int b );
  Rectangle( const int h, const int w );	// left_ = 0, top_ = 0

  void left  ( const int l );
  void top   ( const int t );
  void right ( const int r );
  void bottom( const int b );
  void height( const int h );
  void width ( const int w );
  void add_point( const int row, const int col )
    { if( row > bottom_ ) bottom_ = row; else if( row < top_ ) top_ = row;
      if( col > right_ ) right_ = col;   else if( col < left_ ) left_ = col; }
  void add_rectangle( const Rectangle & re )
    { if( re.left_ < left_ )     left_ = re.left_;
      if( re.top_ < top_ )       top_ = re.top_;
      if( re.right_ > right_ )   right_ = re.right_;
      if( re.bottom_ > bottom_ ) bottom_ = re.bottom_; }

  int left()    const { return left_;   }
  int top()     const { return top_;    }
  int right()   const { return right_;  }
  int bottom()  const { return bottom_; }
  int height()  const { return bottom_ - top_ + 1; }
  int width()   const { return right_ - left_ + 1; }
  int size()    const { return height() * width(); }
  int hcenter() const { return ( left_ + right_ ) / 2; }
  int vcenter() const { return ( top_ + bottom_ ) / 2; }
  int hpos( const int p ) const
    { return left_ + ( ( ( right_ - left_ ) * p ) / 100 ); }
  int vpos( const int p ) const
    { return top_ + ( ( ( bottom_ - top_ ) * p ) / 100 ); }

  bool operator==( const Rectangle & re ) const
    { return ( left_ == re.left_ && top_ == re.top_ && right_ == re.right_ && bottom_ == re.bottom_ ); }
  bool operator!=( const Rectangle & re ) const { return !( *this == re ); }

  bool includes( const Rectangle & re ) const
    { return left_  <= re.left_  && top_    <= re.top_ &&
             right_ >= re.right_ && bottom_ >= re.bottom_; }
  bool includes( const int row, const int col ) const
    { return left_ <= col && right_ >= col && top_ <= row && bottom_ >= row; }
  bool strictly_includes( const Rectangle & re ) const
    { return left_  < re.left_  && top_    < re.top_ &&
             right_ > re.right_ && bottom_ > re.bottom_; }
  bool strictly_includes( const int row, const int col ) const
    { return left_ < col && right_ > col && top_ < row && bottom_ > row; }
  bool includes_hcenter( const Rectangle & re ) const
    { const int hc = re.hcenter(); return left_ <= hc && right_ >= hc; }
  bool includes_vcenter( const Rectangle & re ) const
    { const int vc = re.vcenter(); return top_ <= vc && bottom_ >= vc; }
  bool h_includes( const Rectangle & re ) const
    { return left_ <= re.left_ && right_ >= re.right_; }
  bool h_includes( const int col ) const
    { return left_ <= col && right_ >= col; }
  bool v_includes( const Rectangle & re ) const
    { return top_ <= re.top_ && bottom_ >= re.bottom_; }
  bool v_includes( const int row ) const
    { return top_ <= row && bottom_ >= row; }
  bool h_overlaps( const Rectangle & re ) const
    { return left_ <= re.right_ && right_ >= re.left_; }
  bool v_overlaps( const Rectangle & re ) const
    { return top_ <= re.bottom_ && bottom_ >= re.top_; }
  int  v_overlap_percent( const Rectangle & re ) const;
  bool is_hcentred_in( const Rectangle & re ) const;
  bool is_vcentred_in( const Rectangle & re ) const;
  bool precedes( const Rectangle & re ) const;
  bool h_precedes( const Rectangle & re ) const
    { return hcenter() < re.hcenter(); }
  bool v_precedes( const Rectangle & re ) const;

  int distance( const Rectangle & re ) const
    { return hypoti( h_distance( re ), v_distance( re ) ); }
  int distance( const int row, const int col ) const
    { return hypoti( h_distance( col ), v_distance( row ) ); }
  int h_distance( const Rectangle & re ) const
    { if( re.right_ <= left_ ) return left_ - re.right_;
      if( re.left_ >= right_ ) return re.left_ - right_;
      return 0; }
  int h_distance( const int col ) const
    { if( col <= left_ ) return left_ - col;
      if( col >= right_ ) return col - right_;
      return 0; }
  int v_distance( const Rectangle & re ) const
    { if( re.bottom_ <= top_ ) return top_ - re.bottom_;
      if( re.top_ >= bottom_ ) return re.top_ - bottom_;
      return 0; }
  int v_distance( const int row ) const
    { if( row <= top_ ) return top_ - row;
      if( row >= bottom_ ) return row - bottom_;
      return 0; }

  static int hypoti( const int c1, const int c2 );
  };
