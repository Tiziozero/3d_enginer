package main
import "core:fmt"
import "core:math"
import "vendor:raylib"

SCALE :: 100
WIDTH :: 400
HEIGHT :: 250
log :: proc(a: ..any){}//fmt.print
logln :: proc(a: ..any){}//fmt.println


// get coords to draw on screen based on center 
get_draw_coords :: proc(p :Vec2) -> Vec2 {
    return Vec2{WIDTH / 2 + p.x* SCALE, HEIGHT/2+p.y*SCALE}
}

Vec2 :: raylib.Vector2
sin::math.sin
cos::math.cos

normalize2 :: proc(v: Vec2) -> Vec2 {
    mag := math.sqrt(v.x*v.x + v.y*v.y);
    if mag == 0.0 { return Vec2{0, 0}; }
    return Vec2{v.x / mag, v.y / mag};
}

normalize3 :: proc(v: Vec3) -> Vec3 {
    mag := math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if mag == 0.0 { return Vec3{0, 0, 0}; }
    return Vec3{v.x / mag, v.y / mag, v.z / mag};
}
rotate_point :: proc(p, origin: Vec2, angle_rad: f32) -> Vec2 {
    dx := p.x - origin.x;
    dy := p.y - origin.y;

    cos_a := cos(angle_rad);
    sin_a := sin(angle_rad);

    return Vec2{
        origin.x + dx * cos_a - dy * sin_a,
        origin.y + dx * sin_a + dy * cos_a,
    };
}

rotate_rectangle :: proc(points: [4]Vec2, origin: Vec2, angle_rad: f32) -> [4]Vec2 {
    result: [4]Vec2;
    for i in 0..<4 {
        result[i] = rotate_point(points[i], origin, angle_rad);
    }
    return result;
}

Rect::raylib.Rectangle
Vec3::raylib.Vector3
Color::raylib.Color
Ray::raylib.Ray
Body :: struct {
    body: Rect,
    position: Vec2,
    rotation: f64,
    mass:   f64,
    color: Color,

    velocity: f64,
    angular_velocity:f64,

    acceleration:f64,
    angular_acceleration:f64,

    force: Vec2,
    torque: f64,

    inertia: f64,

    // coefficion of restetution
    cor: f64,
    // coefficion of friction
    cof: f64,

    mesh: []Triangle,
}
dot :: proc(a, b: Vec3) -> f32 {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
cross :: proc(a, b: Vec3) -> Vec3 {
    return Vec3{
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x,
    };
}

Triangle :: [3]Vec3


MollerTrumbore :: proc(v1, v2, v0 :Vec3, ray:Ray) -> (bool, f32) {
    e1 := v1 - v0;
    e2 := v2 - v0;

    p := cross (ray.direction, e2);
    det := dot(p, e1);
    epsilon : f32 = 1e-2;
    if (abs(det) <= epsilon) {
        return false, 0;
    }
    inv_det := 1.0/det;


    t_vec := ray.position - v0;

    u := (dot(t_vec,p)) * inv_det;
    if u < 0|| u > 1 {
        return false,0;
    }

    q := cross(t_vec,e1);

    v := dot(ray.direction, q)*inv_det;
    if v < 0|| v > 1 {
        return false,0;
    }
    if u + v > 1 {
        return false, 0;
    }

    t := dot(e2, q) * inv_det;
    if t > 0 {
        return true, t;
    }

    return false, 0;
}


draw_body :: proc(b: Body) {
    raylib.DrawRectanglePro(
        b.body, b.position, f32(b.rotation), b.color
    )
}

a : quaternion256

WORLD_UP :: Vec3{0,1,0};
Camera ::struct {
    pos: Vec3,
    c_forward: Vec3,
    focal_length: f32,
    c_right: Vec3,
    c_up: Vec3,
    s_center: Vec3, // c_forward * focal_length
    s_width: f32,
    s_height: f32,
    s_horizontal: Vec3, // c_right * s_width
    s_vertical: Vec3, // c_up * s_height

    c_angle_h: f32,
    c_angle_v: f32,
}
update_camera:: proc(c: ^Camera)-> Camera {
    s_ah := sin(c.c_angle_h)
    c_ah := cos(c.c_angle_h)
    s_av := sin(c.c_angle_v)
    c_av := cos(c.c_angle_v)
    d_ := c_av
    c.c_forward = normalize3(Vec3{
        -s_ah * c_av,
        s_av,
        c_ah * c_av,
    });

    c.s_center = c.c_forward * c.focal_length;
    c.c_right = normalize3(cross(WORLD_UP, c.c_forward));
    c.c_up = normalize3(cross(c.c_right,c.c_forward));

    c.s_horizontal = c.c_right * c.s_width
    c.s_vertical = c.c_up * c.s_height
    return c^
}
camera_bottom_left :: proc(c: ^Camera) -> Vec3 {
    return c.s_center - c.s_horizontal * 0.5 - c.s_vertical * 0.5;
}

camera_bottom_right :: proc(c: ^Camera) -> Vec3 {
    return c.s_center + c.s_horizontal * 0.5 - c.s_vertical * 0.5;
}

camera_top_left :: proc(c: ^Camera) -> Vec3 {
    return c.s_center - c.s_horizontal * 0.5 + c.s_vertical * 0.5;
}

camera_top_right :: proc(c: ^Camera) -> Vec3 {
    return c.s_center + c.s_horizontal * 0.5 + c.s_vertical * 0.5;
}


main :: proc() {
    logln("Hello, Oding")
    raylib.InitWindow(WIDTH, HEIGHT, "Engine 1")

    c:= Camera {
        Vec3{0,0,0},
        Vec3{0,0,1},
        2,
        Vec3{0,0,0},
        Vec3{0,0,0},
        Vec3{0,0,0},
        WIDTH/SCALE,
        HEIGHT/SCALE,
        Vec3{0,0,0},
        Vec3{0,0,0},
        0,
        0,
    };
    c =update_camera(&c);
    logln(c)
    bodies := []Body{
        Body{
            mesh= []Triangle{
                // CUBE 1 - Center cube (your original)
                // Front face (z = +0.25)
                {Vec3{-0.25, -0.25, 0.25}, Vec3{0.25, -0.25, 0.25}, Vec3{0.25, 0.25, 0.25}},
                {Vec3{-0.25, -0.25, 0.25}, Vec3{0.25, 0.25, 0.25}, Vec3{-0.25, 0.25, 0.25}},
                // Back face (z = -0.25)
                {Vec3{0.25, -0.25, -0.25}, Vec3{-0.25, -0.25, -0.25}, Vec3{-0.25, 0.25, -0.25}},
                {Vec3{0.25, -0.25, -0.25}, Vec3{-0.25, 0.25, -0.25}, Vec3{0.25, 0.25, -0.25}},
                // Left face (x = -0.25)
                {Vec3{-0.25, -0.25, -0.25}, Vec3{-0.25, -0.25, 0.25}, Vec3{-0.25, 0.25, 0.25}},
                {Vec3{-0.25, -0.25, -0.25}, Vec3{-0.25, 0.25, 0.25}, Vec3{-0.25, 0.25, -0.25}},
                // Right face (x = +0.25)
                {Vec3{0.25, -0.25, 0.25}, Vec3{0.25, -0.25, -0.25}, Vec3{0.25, 0.25, -0.25}},
                {Vec3{0.25, -0.25, 0.25}, Vec3{0.25, 0.25, -0.25}, Vec3{0.25, 0.25, 0.25}},
                // Top face (y = +0.25)
                {Vec3{-0.25, 0.25, 0.25}, Vec3{0.25, 0.25, 0.25}, Vec3{0.25, 0.25, -0.25}},
                {Vec3{-0.25, 0.25, 0.25}, Vec3{0.25, 0.25, -0.25}, Vec3{-0.25, 0.25, -0.25}},
                // Bottom face (y = -0.25)
                {Vec3{-0.25, -0.25, -0.25}, Vec3{0.25, -0.25, -0.25}, Vec3{0.25, -0.25, 0.25}},
                {Vec3{-0.25, -0.25, -0.25}, Vec3{0.25, -0.25, 0.25}, Vec3{-0.25, -0.25, 0.25}},
            },
            color=raylib.RED,
        },
        Body{
            mesh= []Triangle{

                // CUBE 2 - Left cube (offset by -1.0 in X)
                // Front face
                {Vec3{-1.25, -0.25, 0.25}, Vec3{-0.75, -0.25, 0.25}, Vec3{-0.75, 0.25, 0.25}},
                {Vec3{-1.25, -0.25, 0.25}, Vec3{-0.75, 0.25, 0.25}, Vec3{-1.25, 0.25, 0.25}},
                // Back face
                {Vec3{-0.75, -0.25, -0.25}, Vec3{-1.25, -0.25, -0.25}, Vec3{-1.25, 0.25, -0.25}},
                {Vec3{-0.75, -0.25, -0.25}, Vec3{-1.25, 0.25, -0.25}, Vec3{-0.75, 0.25, -0.25}},
                // Left face
                {Vec3{-1.25, -0.25, -0.25}, Vec3{-1.25, -0.25, 0.25}, Vec3{-1.25, 0.25, 0.25}},
                {Vec3{-1.25, -0.25, -0.25}, Vec3{-1.25, 0.25, 0.25}, Vec3{-1.25, 0.25, -0.25}},
                // Right face
                {Vec3{-0.75, -0.25, 0.25}, Vec3{-0.75, -0.25, -0.25}, Vec3{-0.75, 0.25, -0.25}},
                {Vec3{-0.75, -0.25, 0.25}, Vec3{-0.75, 0.25, -0.25}, Vec3{-0.75, 0.25, 0.25}},
                // Top face
                {Vec3{-1.25, 0.25, 0.25}, Vec3{-0.75, 0.25, 0.25}, Vec3{-0.75, 0.25, -0.25}},
                {Vec3{-1.25, 0.25, 0.25}, Vec3{-0.75, 0.25, -0.25}, Vec3{-1.25, 0.25, -0.25}},
                // Bottom face
                {Vec3{-1.25, -0.25, -0.25}, Vec3{-0.75, -0.25, -0.25}, Vec3{-0.75, -0.25, 0.25}},
                {Vec3{-1.25, -0.25, -0.25}, Vec3{-0.75, -0.25, 0.25}, Vec3{-1.25, -0.25, 0.25}},
            },color=raylib.BLUE,
        },
        Body{
            mesh= []Triangle{

                // CUBE 3 - Right cube (offset by +1.0 in X)
                // Front face
                {Vec3{0.75, -0.25, 0.25}, Vec3{1.25, -0.25, 0.25}, Vec3{1.25, 0.25, 0.25}},
                {Vec3{0.75, -0.25, 0.25}, Vec3{1.25, 0.25, 0.25}, Vec3{0.75, 0.25, 0.25}},
                // Back face
                {Vec3{1.25, -0.25, -0.25}, Vec3{0.75, -0.25, -0.25}, Vec3{0.75, 0.25, -0.25}},
                {Vec3{1.25, -0.25, -0.25}, Vec3{0.75, 0.25, -0.25}, Vec3{1.25, 0.25, -0.25}},
                // Left face
                {Vec3{0.75, -0.25, -0.25}, Vec3{0.75, -0.25, 0.25}, Vec3{0.75, 0.25, 0.25}},
                {Vec3{0.75, -0.25, -0.25}, Vec3{0.75, 0.25, 0.25}, Vec3{0.75, 0.25, -0.25}},
                // Right face
                {Vec3{1.25, -0.25, 0.25}, Vec3{1.25, -0.25, -0.25}, Vec3{1.25, 0.25, -0.25}},
                {Vec3{1.25, -0.25, 0.25}, Vec3{1.25, 0.25, -0.25}, Vec3{1.25, 0.25, 0.25}},
                // Top face
                {Vec3{0.75, 0.25, 0.25}, Vec3{1.25, 0.25, 0.25}, Vec3{1.25, 0.25, -0.25}},
                {Vec3{0.75, 0.25, 0.25}, Vec3{1.25, 0.25, -0.25}, Vec3{0.75, 0.25, -0.25}},
                // Bottom face
                {Vec3{0.75, -0.25, -0.25}, Vec3{1.25, -0.25, -0.25}, Vec3{1.25, -0.25, 0.25}},
                {Vec3{0.75, -0.25, -0.25}, Vec3{1.25, -0.25, 0.25}, Vec3{0.75, -0.25, 0.25}},
            },color=raylib.GREEN,
        },
        Body{
            mesh=[]Triangle{

                // CUBE 4 - Top cube (offset by +1.0 in Y)
                // Front face
                {Vec3{-0.25, 0.75, 0.25}, Vec3{0.25, 0.75, 0.25}, Vec3{0.25, 1.25, 0.25}},
                {Vec3{-0.25, 0.75, 0.25}, Vec3{0.25, 1.25, 0.25}, Vec3{-0.25, 1.25, 0.25}},
                // Back face
                {Vec3{0.25, 0.75, -0.25}, Vec3{-0.25, 0.75, -0.25}, Vec3{-0.25, 1.25, -0.25}},
                {Vec3{0.25, 0.75, -0.25}, Vec3{-0.25, 1.25, -0.25}, Vec3{0.25, 1.25, -0.25}},
                // Left face
                {Vec3{-0.25, 0.75, -0.25}, Vec3{-0.25, 0.75, 0.25}, Vec3{-0.25, 1.25, 0.25}},
                {Vec3{-0.25, 0.75, -0.25}, Vec3{-0.25, 1.25, 0.25}, Vec3{-0.25, 1.25, -0.25}},
                // Right face
                {Vec3{0.25, 0.75, 0.25}, Vec3{0.25, 0.75, -0.25}, Vec3{0.25, 1.25, -0.25}},
                {Vec3{0.25, 0.75, 0.25}, Vec3{0.25, 1.25, -0.25}, Vec3{0.25, 1.25, 0.25}},
                // Top face
                {Vec3{-0.25, 1.25, 0.25}, Vec3{0.25, 1.25, 0.25}, Vec3{0.25, 1.25, -0.25}},
                {Vec3{-0.25, 1.25, 0.25}, Vec3{0.25, 1.25, -0.25}, Vec3{-0.25, 1.25, -0.25}},
                // Bottom face
                {Vec3{-0.25, 0.75, -0.25}, Vec3{0.25, 0.75, -0.25}, Vec3{0.25, 0.75, 0.25}},
                {Vec3{-0.25, 0.75, -0.25}, Vec3{0.25, 0.75, 0.25}, Vec3{-0.25, 0.75, 0.25}},
            },color=raylib.GRAY
        },
        Body{
            mesh=[]Triangle{

                // CUBE 5 - Back cube (offset by -1.0 in Z)
                // Front face
                {Vec3{-0.25, -0.25, -0.75}, Vec3{0.25, -0.25, -0.75}, Vec3{0.25, 0.25, -0.75}},
                {Vec3{-0.25, -0.25, -0.75}, Vec3{0.25, 0.25, -0.75}, Vec3{-0.25, 0.25, -0.75}},
                // Back face
                {Vec3{0.25, -0.25, -1.25}, Vec3{-0.25, -0.25, -1.25}, Vec3{-0.25, 0.25, -1.25}},
                {Vec3{0.25, -0.25, -1.25}, Vec3{-0.25, 0.25, -1.25}, Vec3{0.25, 0.25, -1.25}},
                // Left face
                {Vec3{-0.25, -0.25, -1.25}, Vec3{-0.25, -0.25, -0.75}, Vec3{-0.25, 0.25, -0.75}},
                {Vec3{-0.25, -0.25, -1.25}, Vec3{-0.25, 0.25, -0.75}, Vec3{-0.25, 0.25, -1.25}},
                // Right face
                {Vec3{0.25, -0.25, -0.75}, Vec3{0.25, -0.25, -1.25}, Vec3{0.25, 0.25, -1.25}},
                {Vec3{0.25, -0.25, -0.75}, Vec3{0.25, 0.25, -1.25}, Vec3{0.25, 0.25, -0.75}},
                // Top face
                {Vec3{-0.25, 0.25, -0.75}, Vec3{0.25, 0.25, -0.75}, Vec3{0.25, 0.25, -1.25}},
                {Vec3{-0.25, 0.25, -0.75}, Vec3{0.25, 0.25, -1.25}, Vec3{-0.25, 0.25, -1.25}},
                // Bottom face
                {Vec3{-0.25, -0.25, -1.25}, Vec3{0.25, -0.25, -1.25}, Vec3{0.25, -0.25, -0.75}},
                {Vec3{-0.25, -0.25, -1.25}, Vec3{0.25, -0.25, -0.75}, Vec3{-0.25, -0.25, -0.75}},
            },color=raylib.YELLOW
        },
    }
    raylib.DisableCursor()

    k :f32 = 0
    for (!raylib.WindowShouldClose()) {
        dt := raylib.GetFrameTime();
        if raylib.IsKeyDown(raylib.KeyboardKey.LEFT ) do c.pos -= c.c_right*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.RIGHT) do c.pos += c.c_right*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.A) do c.pos -= c.c_right*dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.D) do c.pos += c.c_right*dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.UP) do c.pos += c.c_forward * dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.DOWN) do c.pos -= c.c_forward * dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.W) do c.pos += c.c_forward * dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.S) do c.pos -= c.c_forward * dt;

        if raylib.IsKeyDown(raylib.KeyboardKey.L) do c.c_angle_h += 0.01;
        if raylib.IsKeyDown(raylib.KeyboardKey.R) do c.c_angle_h -= 0.01;
        if raylib.IsKeyDown(raylib.KeyboardKey.H) do c.c_angle_h -= 0.01;

        if raylib.IsKeyDown(raylib.KeyboardKey.SPACE) do c.pos.y -=  dt;
        if raylib.IsKeyDown(raylib.KeyboardKey.LEFT_SHIFT) do c.pos.y +=  dt;

        mouse_delta := raylib.GetMouseDelta()
        if mouse_delta.x != 0 do c.c_angle_h -= mouse_delta.x/WIDTH * 2*3.1415693;
        if mouse_delta.y != 0 do c.c_angle_v -= mouse_delta.y/HEIGHT * 2*3.1415693;


        update_camera(&c);
        tl := camera_top_left(&c);
        log(c.pos);
        log(c.c_forward);
        log(c.c_angle_h);
        log(tl);
        logln();
        raylib.BeginDrawing();
        raylib.ClearBackground(raylib.BLACK);

        for body in bodies {
            mesh := body.mesh
            render_mesh2(&c, mesh, body.color);
            draw_mesh(mesh)
        }
        draw_camera(c)
        draw_screen_calcs(c);
        raylib.EndDrawing();
    }
    raylib.EnableCursor()
}


draw_camera :: proc(c:Camera) {
    scale :f32= SCALE/2
    x := i32(c.pos.x*scale+200)
    y := i32(c.pos.z*scale+200)

    raylib.DrawCircle(x, y, 5, raylib.BLUE);
    dx :i32= x + i32(c.c_forward.x*f32(scale))
    dy :i32= (y + i32(c.c_forward.z*f32(scale)))
    raylib.DrawLine(x,y,dx,dy,raylib.RED);
}

draw_mesh :: proc(mesh: []Triangle) {
    scale :f32= SCALE/2
    for t in mesh {
        for v in t {
            x := i32(v.x*scale+200)
            y := i32(v.z*scale+200)
            raylib.DrawCircle(x, y, 1, raylib.RED);
        }
    }

}
render_mesh :: proc(c:^Camera, mesh: []Triangle) {
    tl := camera_top_left(c)
    k :f32= 1/SCALE;
    fin := false
    for i :i32= 0; i < WIDTH; i+= 1 {
        for j :i32= 0; j < HEIGHT; j+=1 {
            r := Ray{
                c.pos,
                tl + Vec3{f32(i)/SCALE, f32(j)/SCALE,0},
            }

            closest_t := math.inf_f32(1)
            hit_something := false
            for tr in mesh {
                if hit, t := MollerTrumbore(tr[0],tr[1],tr[2], r); hit && t < closest_t {
                    closest_t = t
                    hit_something = true
                    raylib.DrawPixel(i,j,raylib.GRAY);
                }
            }
        }
    }
}
render_mesh2 :: proc(c:^Camera, mesh: []Triangle, color: Color) {
    tl := camera_top_left(c)

    for i :i32= 0; i < WIDTH; i+= 1 {
        for j :i32= 0; j < HEIGHT; j+=1 {
            // Calculate the target point on the screen plane
            u := f32(i) / f32(WIDTH)   // 0 to 1
            v := f32(j) / f32(HEIGHT)  // 0 to 1
            target := tl + c.s_horizontal * u + c.s_vertical * v;
            direction := normalize3(target - c.pos)
            r := Ray{
                c.pos,
                direction,
            }

            closest_t := math.inf_f32(1)
            hit_something := false
            for tr in mesh {
                if hit, t := MollerTrumbore(tr[0],tr[1],tr[2], r); hit && t < closest_t && t > 0 {
                    closest_t = t
                    hit_something = true
                }
            }

            if hit_something {
                raylib.DrawPixel(i, j, color)
            }
        }
    }
}


draw_screen_calcs  :: proc(c:Camera) {
    cam := c
    scale := SCALE/3
    x0 :i32= WIDTH / 2
    y0 :i32= HEIGHT / 2
    m := []Vec3{
        camera_top_left(&cam),
        camera_top_right(&cam),
        camera_bottom_left(&cam),
        camera_bottom_right(&cam),
    }
    for i in m {
        log(i)
        raylib.DrawCircle(x0+i32(i.x*SCALE), y0+i32(i.y*SCALE), 2,raylib.WHITE);
    }
    logln()
    raylib.DrawCircle(x0+i32(c.s_center.x*SCALE), y0+i32(c.s_center.y*SCALE), 2,raylib.WHITE);

}
