#========================
Date conversion
=========================#

function datetime_to_days(datetime::Dates.DateTime)
	# Converts datetime object to days since 01.01.1970 as Float64.
	(datetime-Dates.DateTime(1970)).value/1.0e3/3600.0/24.0
end
precompile(datetime_to_days, (Dates.DateTime,))

function days_to_datetime(days)
	#Converts days since 01.01.1970 to DateTime object.
	Dates.DateTime(1970)+Dates.Nanosecond(floor(days*24.0*3600.0*1.0e9))
end
precompile(days_to_datetime, (Float64,))
precompile(days_to_datetime, (Int64,))