from django.shortcuts import render
from django import forms
from .tasks import (
    sectionb,
    hohmann_transfer,
    jd_lst,
    sectiona,
    ref_frame_and_coords_transfer,
    sectionf
)

class SectionB1Form(forms.Form):
    r0 = forms.CharField(
        label="R0", widget=forms.TextInput(attrs={"placeholder": "1 2 3"})
    )
    v0 = forms.CharField(
        label="V0", widget=forms.TextInput(attrs={"placeholder": "1 2 3"})
    )
    t = forms.IntegerField(label="t")


class SectionB2Form(forms.Form):
    r = forms.CharField(
        label="R", widget=forms.TextInput(attrs={"placeholder": "1 2 3"})
    )
    v = forms.CharField(
        label="V", widget=forms.TextInput(attrs={"placeholder": "1 2 3"})
    )


class SectionB3Form(forms.Form):
    h = forms.FloatField(label="H")
    e = forms.FloatField(label="E")
    ra = forms.FloatField(label="RA")
    incl = forms.FloatField(label="Incl")
    w = forms.FloatField(label="W")
    ta = forms.FloatField(label="Ta")


class Hoffman(forms.Form):
    ri = forms.FloatField(label="Ri")
    rf = forms.FloatField(label="Rf")
    mu = forms.FloatField(label="Mu")


class Julian(forms.Form):
    year = forms.FloatField(label=" Year")
    month = forms.FloatField(label="Month")
    day = forms.FloatField(label="Day")


class LST(forms.Form):
    jul_date = forms.FloatField(label="Julian Date")
    longitude = forms.FloatField(label="Longitude")
    utc_time = forms.FloatField(label="UTC_Time")
class SectionF(forms.Form):
    x0=  forms.FloatField(label="X0(Initial X component)")
    y0=  forms.FloatField(label="Y0(Initial Y component)")
    vx0 = forms.FloatField(label="Vx0(Initial X velocity component)")
    vy0 = forms.FloatField(label="Vy0(Initial y velocity component)")
    t0 = forms.FloatField(label="T0(Initial Time)")
    tmax = forms.FloatField(label="Tmax(Final Time)")
    dt = forms.FloatField(label="Step Size")


class Kepler_hyp(forms.Form):
    m = forms.FloatField(label="M")
    e = forms.FloatField(label="E")


class Kepler_Eqn(forms.Form):
    m = forms.FloatField(label="M")
    e = forms.FloatField(label="E")


class SectionC1Form(forms.Form):
    inertial_coords = forms.CharField(
        label="Inertial Coordinates",
        widget=forms.TextInput(attrs={"placeholder": "1 2 3"}),
    )
    rotation_angle = forms.FloatField(label="Rotation Angle")


class SectionC2Form(forms.Form):
    rotating_coords = forms.CharField(
        label="Rotating Coordinates",
        widget=forms.TextInput(attrs={"placeholder": "1 2 3"}),
    )
    rotation_angle = forms.FloatField(label="Rotation Angle")


# Create your views here.
def home(request):
    return render(request, "home.html")


def section_b_1(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionB1Form(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            t = int(form.cleaned_data["t"])
            r0 = list(map(float, form.cleaned_data["r0"].split(" ")))
            v0 = list(map(float, form.cleaned_data["v0"].split(" ")))
            r, v = sectionb.main1(r0, v0, t)
            context["r_val"] = ", ".join(map(str, r))
            context["v_val"] = ", ".join(map(str, v))

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionB1Form()

        context["r_val"] = ""
        context["v_val"] = ""

    context["form"] = form

    return render(request, "section_b_1.html", context)


def section_b_2(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionB2Form(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            r = list(map(float, form.cleaned_data["r"].split(" ")))
            v = list(map(float, form.cleaned_data["v"].split(" ")))
            h, incl, ra, e, w, ta, a, t = sectionb.main2(r, v)
            context["h"] = h
            context["incl"] = incl
            context["ra"] = ra
            context["e"] = e
            context["w"] = w
            context["ta"] = ta
            context["a"] = a
            context["t"] = t

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionB2Form()

        context["h"] = ""
        context["incl"] = ""
        context["ra"] = ""
        context["e"] = ""
        context["w"] = ""
        context["ta"] = ""
        context["a"] = ""
        context["t"] = ""

    context["form"] = form

    return render(request, "section_b_2.html", context)


def section_b_3(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionB3Form(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            h = float(form.cleaned_data["h"])
            e = float(form.cleaned_data["e"])
            ra = float(form.cleaned_data["ra"])
            incl = float(form.cleaned_data["incl"])
            w = float(form.cleaned_data["w"])
            ta = float(form.cleaned_data["ta"])
            r, v = sectionb.main3(h, e, ra, incl, w, ta)
            context["r_val"] = ", ".join(map(str, r))
            context["v_val"] = ", ".join(map(str, v))

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionB3Form()

        context["r_val"] = ""
        context["v_val"] = ""

    context["form"] = form

    return render(request, "section_b_3.html", context)


def hoffman(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = Hoffman(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            ri = float(form.cleaned_data["ri"])
            rf = float(form.cleaned_data["rf"])
            mu = float(form.cleaned_data["mu"])

            ans = hohmann_transfer.hohmann_transfer(ri, rf, mu)
            context["semi_major"] = ans.semi_major_axis
            context["v_peri"] = ans.v_periapsis
            context["v_apo"] = ans.v_apoapsis
            context["vi"] = ans.vi
            context["vf"] = ans.vf
            context["del_vi"] = ans.delta_vi
            context["del_vf"] = ans.delta_vf
            context["tp"] = ans.time_period

    # if a GET (or any other method) we'll create a blank form
    else:
        form = Hoffman()
        context["semi_major"] = ""
        context["v_peri"] = ""
        context["v_apo"] = ""
        context["vi"] = ""
        context["vf"] = ""
        context["del_vi"] = ""
        context["del_vf"] = ""
        context["tp"] = ""
    context["form"] = form

    return render(request, "hoffman.html", context)


def julian(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = Julian(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            year = float(form.cleaned_data["year"])
            month = float(form.cleaned_data["month"])
            day = float(form.cleaned_data["day"])

            ans = jd_lst.calculate_julian_date(year, month, day)
            context["jd"] = ans

    # if a GET (or any other method) we'll create a blank form
    else:
        form = Julian()

        context["jd"] = ""

    context["form"] = form

    return render(request, "julian.html", context)


def lst(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = LST(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            jul_date = float(form.cleaned_data["jul_date"])
            longitude = float(form.cleaned_data["longitude"])
            utc_time = float(form.cleaned_data["utc_time"])

            ans = jd_lst.calculate_lst(jul_date, longitude, utc_time)
            context["lst"] = ans

    # if a GET (or any other method) we'll create a blank form
    else:
        form = LST()

        context["lst"] = ""

    context["form"] = form

    return render(request, "lst.html", context)

def section_a(request):
    return render(request, "section_a.html")


def kepler_hyp(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = Kepler_hyp(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            m = float(form.cleaned_data["m"])
            e = float(form.cleaned_data["e"])
            tol = 1e-9
            ans = sectiona.kepler_hyperbolic(e, m, tol, 1000)
            context["e_m"] = ans

    # if a GET (or any other method) we'll create a blank form
    else:
        form = Kepler_hyp()

        context["e_m"] = ""

    context["form"] = form

    return render(request, "kepler_hyp.html", context)


def kepler_eqn(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = Kepler_Eqn(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            m = float(form.cleaned_data["m"])
            e = float(form.cleaned_data["e"])
            tol = 1e-9
            ans = sectiona.kepler_equation(e, m, tol, 1000)
            context["em"] = ans

    # if a GET (or any other method) we'll create a blank form
    else:
        form = Kepler_Eqn()

        context["em"] = ""

    context["form"] = form

    return render(request, "kepler_eqn.html", context)


def section_c_1(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionC1Form(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            inertial_coords = list(
                map(float, form.cleaned_data["inertial_coords"].split(" "))
            )
            rotation_angle = float(form.cleaned_data["rotation_angle"])
            ans = ref_frame_and_coords_transfer.inertial_to_rotating_frame(
                inertial_coords, rotation_angle
            )
            context["rotation_coords"] = ", ".join(map(str, ans))

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionC1Form()

        context["rotation_coords"] = ""

    context["form"] = form

    return render(request, "section_c_1.html", context)


def section_c_2(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionC2Form(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL
            rotating_coords = list(
                map(float, form.cleaned_data["rotating_coords"].split(" "))
            )
            rotation_angle = float(form.cleaned_data["rotation_angle"])
            ans = ref_frame_and_coords_transfer.rotating_to_inertial_frame(
                rotating_coords, rotation_angle
            )
            context["inertial_coords"] = ", ".join(map(str, ans))

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionC2Form()

        context["inertial_coords"] = ""

    context["form"] = form

    return render(request, "section_c_2.html", context)


def section_f(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionF(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            x0= float(form.cleaned_data["x0"])
            y0 = float(form.cleaned_data["y0"])
            vx0 = float(form.cleaned_data["vx0"])
            vy0 = float(form.cleaned_data["vy0"])
            t0 = float(form.cleaned_data["t0"])
            tmax = float(form.cleaned_data["tmax"])
            dt = float(form.cleaned_data["dt"])
            ans = sectionf.calc(x0,y0,vx0,vy0,t0,tmax,dt)
            context["traj"] = ",\n".join(map(str, ans))

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionF()

        context["traj"] = ""


    context["form"] = form

    return render(request, "section_f.html", context)
